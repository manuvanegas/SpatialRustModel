export step_model!

function step_model!(model::SpatialRustABM)
    pre_step!(model)
    grow_shades!(model.current, model.mngpars.shade_g_rate, model.mngpars.max_shade)
    coffee_step!(model)
    rust_step!(model)
    farmer_step!(model)
    return nothing
end

## "Step" functions

function pre_step!(model::SpatialRustABM)
    # update day counters
    model.current.days += 1
    t = model.current.ticks += 1

    # update weather conditions from Weather data
    model.current.rain = model.weather.rain_data[t]
    @inbounds model.current.wind = model.weather.wind_data[t]
    @inbounds model.current.temperature = model.weather.temp_data[t]

    # spore outpour decay, then outpour can return spores to the farm if windy
    model.outpour .*= 0.9
    if model.current.wind
        model.current.wind_h = rand(model.rng) * 360.0
    end

    # used in GA, in case the evolved strategy was too effective and all Rusts were eliminated:
    # if t % 30 == 0 && sum(activeRust, model.agents) == 0 && sum(active, model.agents) > 400
    #     reintroduce_rusts!(model, 5)
    # end
    return nothing
end

function coffee_step!(model::SpatialRustABM)
    prod_cycle_d = mod1(model.current.days, 365)
    pars = model.coffeepars

    if pars.veg_d < pars.rep_d
        vegd = pars.veg_d
        repd = pars.rep_d
        cycled = prod_cycle_d
    else
        vegd = - pars.veg_d
        repd = - pars.rep_d
        cycled = - prod_cycle_d
    end

    if cycled != repd
        if vegd <= cycled < repd
            growth = veg_growth!
        else
            growth = rep_growth!
        end
        for cof in model.agents
            if cof.exh_countdown == 0
                sl = update_sunlight!(cof, model.shade_map, model.current.ind_shade)
                growth(cof, pars, sl)
            elseif cof.exh_countdown > 1
                cof.exh_countdown -= 1
            else
                sl = update_sunlight!(cof, model.shade_map, model.current.ind_shade)
                regrow!(cof, sl, model.farm_map)
            end
        end
    else
        commit_dist = Normal(pars.res_commit, 0.01)
        for cof in model.agents
            if cof.exh_countdown == 0
                sl = update_sunlight!(cof, model.shade_map, model.current.ind_shade)
                veg_growth!(cof, pars, sl)
                cof.production = max(0.0, rand(model.rng, commit_dist) * cof.sunlight * cof.veg * cof.storage)
            elseif cof.exh_countdown > 1
                cof.exh_countdown -= 1
            else
                sl = update_sunlight!(cof, model.shade_map, model.current.ind_shade)
                regrow!(cof, sl, model.farm_map)
            end
        end
    end
    return nothing
end

function rust_step!(model::SpatialRustABM)
    # three independent conditions: fungicide in effect? rainy? windy? -> 8 options
    # not the most clean/maintainable implementation, but I was prioritizing min sim time
    if model.current.fungicide == 0
        if model.current.rain
            if model.current.wind
                rust_step_schedule(model, 1.0, 0, 1.0, r_germinate!, grow_rust!, disperse_rain!, disperse_wind!)
                outside_spores!(model)
            else
                rust_step_schedule(model, 1.0, 0, 1.0, r_germinate!, grow_rust!, disperse_rain!, dummy_disp)
            end
        else
            dry_spo = model.rustpars.pdry_spo
            if model.current.wind
                rust_step_schedule(model, 1.0, 0, dry_spo, nr_germinate!, grow_rust!, dummy_disp, disperse_wind!)
                outside_spores!(model)
            else
                rust_step_schedule(model, 1.0, 0, dry_spo, nr_germinate!, grow_rust!, dummy_disp, dummy_disp)
            end
        end
    else
        f_day = model.current.fungicide
        fung_inf = model.rustpars.fung_inf
        if model.current.rain
            if model.current.wind
                rust_step_schedule(model, fung_inf, f_day, 1.0, r_germinate!, grow_f_rust!, disperse_rain!, disperse_wind!)
                outside_spores!(model)
            else
                rust_step_schedule(model, fung_inf, f_day, 1.0, r_germinate!, grow_f_rust!, disperse_rain!, dummy_disp)
            end
        else
            dry_spo = model.rustpars.pdry_spo
            if model.current.wind
                rust_step_schedule(model, fung_inf, f_day, dry_spo, nr_germinate!, grow_f_rust!, dummy_disp, disperse_wind!)
                outside_spores!(model)
            else
                rust_step_schedule(model, fung_inf, f_day, dry_spo, nr_germinate!, grow_f_rust!, dummy_disp, dummy_disp)
            end
        end
    end
    return nothing
end

function rust_step_schedule(model::SpatialRustABM, f_inf::Float64, f_day::Int, rain_spo::Float64, germinate_f::Function, grow_f::Function,
    rain_dispersal::Function, wind_dispersal::Function)
    
    rusts = Iterators.filter(r -> r.rusted, model.agents)

    for rust in rusts
        if any(rust.spores)
            spore_area = sum(last(p) for p in pairs(rust.areas) if rust.spores[first(p)]) * model.rustpars.spore_pct * (1.0 + rust.sunlight)
            rain_dispersal(model, rust, spore_area)
            wind_dispersal(model, rust, spore_area)
        end
        
        local_temp = model.current.temperature - (model.rustpars.temp_cooling * (1.0 - rust.sunlight))

        if rust.n_lesions > 0
            grow_f(rust, model.rng, model.rustpars, local_temp, rain_spo, f_day)
        
            if losttrack(rust.areas) || !isfinite(rust.storage)
                model.current.withinbounds = false
                break
            end

            germinate_f(rust, model.rng, model.rustpars, local_temp, f_inf)

            parasitize!(rust, model.rustpars, model.farm_map)
        else
            germinate_f(rust, model.rng, model.rustpars, local_temp, f_inf)
        end
    end

    # Update happens in a second loop because first all rusts have had to (try to) disperse
    for rust in rusts
        update_rust!(rust, model.rustpars.viab_loss)
    end
    return nothing
end

function losttrack(as)
    any(a -> (!isfinite(a) || a < -0.1 || a > 25.9), as)
end

function farmer_step!(model)
    doy = mod1(model.current.days, 365)

    if doy == model.mngpars.harvest_day
        harvest!(model)
    end

    prune_i = findfirst(==(doy), model.mngpars.prune_sch)
    if !isnothing(prune_i)
        prune_shades!(model, model.mngpars.post_prune[prune_i])
    end

    if model.current.days % model.mngpars.inspect_period == 0
        inspect!(model)
    end
    
    if model.current.fungicide > 0
        if model.current.fungicide > 30
            model.current.fungicide = 0
        else
            model.current.fungicide += 1
        end
    elseif doy != model.mngpars.harvest_day && model.current.fung_count < model.mngpars.max_fung_sprayings
        if doy in model.mngpars.fungicide_sch
            fungicide!(model)
        elseif model.current.obs_incidence > model.mngpars.incidence_thresh
            fs = model.mngpars.fung_stratg
            if fs == :incd || fs == :cal_incd || (fs == :flor && doy > @inbounds model.mngpars.fungicide_sch[2])
                fungicide!(model)
            end
        end
    end
    
    model.current.shadeacc += model.current.ind_shade
    
    return nothing
end
