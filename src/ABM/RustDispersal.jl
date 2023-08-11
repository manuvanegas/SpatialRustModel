# Spore dispersal and deposition

function disperse_rain!(model::SpatialRustABM, rust::Coffee, spores::Float64)
    map = model.farm_map
    pars = model.rustpars
    pos = rust.pos
    d_mod = (4.0 - 4.0 * pars.diff_splash) * (rust.sunlight - 0.5)^2.0 + pars.diff_splash
    exp_dist = Exponential(pars.rain_dst)
    splashed = rand(model.rng, Poisson(spores))
    for _ in 1:splashed
        distance = rand(model.rng, exp_dist) * d_mod
        if distance < 1.0 #self-infected
            rust.newdeps += 1.0
        else
            # follow splash and return: Tuple > 0 -> Coffee pos, < 0 -> outpour direction (see setup for mappings), 0 -> nothing
            fin_pos = splash(model.rng, pos, distance, rand(model.rng) * 360.0, map, pars)
            if any(fin_pos .> 0) && 
                (c = (@inbounds model[id_in_position(fin_pos,model)])).exh_countdown == 0
                c.newdeps += 1.0
                c.rusted = true
            elseif any(fin_pos .< 0) 
                model.outpour[sum(fin_pos .* (-3,-1))] += 1.0
            end
        end
    end
end

function disperse_wind!(model::SpatialRustABM, rust::Coffee, spores::Float64)
    shading = @inbounds model.shade_map[rust.pos...]
    rustpars = model.rustpars
    w_distance = rand(model.rng, Exponential(rustpars.wind_dst)) * (1 + rust.sunlight * rustpars.diff_wind)
    lifted = rand(model.rng, Poisson(spores * (1.0 - 0.5 * shading)))
    if w_distance < 1.0
        for _ in 1:lifted
            rust.newdeps += 1.0
        end
    else
        wind_h = model.current.wind_h - 15.0
        for _ in 1:lifted
            fin_pos = gust(model.rng, rust.pos, w_distance,
            (wind_h + (rand(model.rng) * 30.0)),
            model.farm_map, model.shade_map, rustpars)
            if any(fin_pos .> 0) && 
                (c = (@inbounds model[id_in_position(fin_pos,model)])).exh_countdown == 0
                c.newdeps += 1.0
                c.rusted = true
            elseif any(fin_pos .< 0) 
                model.outpour[sum(fin_pos .* (-3,-1))] += 1.0
            end
        end
    end
end

dummy_disp(model::SpatialRustABM, rust::Coffee, spores::Float64) = nothing

function splash(rng, pos::NTuple{2,Int}, dist::Float64, heading::Float64, farm_map::Array{Int, 2}, rustpars::RustPars)
    ca = cosd(heading)
    co = sind(heading)
    side = rustpars.map_side
    prob_block = rustpars.tree_block
    px, py = pos

    traveled = 1.0
    onx = 0
    ony = 0
    advanced = false

    while traveled <= dist
        newx = px + floor(Int, ca * traveled)
        newy = py + floor(Int, co * traveled)
        if newx != onx
            onx = newx
            advanced = true
        end
        if newy != ony
            ony = newy
            advanced = true
        end
        if advanced
            withinboundsx = ifelse(newx < 1, -1, ifelse(newx > side, -2, 0))
            withinboundsy = ifelse(newy < 1, -1, ifelse(newy > side, -2, 0))
            if withinboundsx < 0 || withinboundsy < 0
                return (withinboundsx, withinboundsy)
            else
                if @inbounds (id = farm_map[newx, newy]) == 1 && (rand(rng) < prob_block)
                    return (newx, newy)
                elseif id == 2 && (rand(rng) < prob_block)
                    return (0,0)
                end
            end
        end
        advanced = false
        traveled += 0.5
    end

    if @inbounds farm_map[onx, ony] == 1
        (onx, ony)
    else
        return (0,0)
    end
end

function gust(rng, pos::NTuple{2,Int}, dist::Float64, heading::Float64, farm_map::Array{Int,2}, shade_map::Array{Float64, 2}, rustpars::RustPars)
    ca = cosd(heading)
    co = sind(heading)
    side = rustpars.map_side
    prob_block = rustpars.shade_block
    px, py = pos

    notblocked = true
    traveled = 1.0
    onx = 0
    ony = 0
    advanced = false

    while traveled <= dist
        newx = px + floor(Int, ca * traveled)
        newy = py + floor(Int, co * traveled)
        if newx != onx
            onx = newx
            advanced = true
        end
        if newy != ony
            ony = newy
            advanced = true
        end
        if advanced
            withinboundsx = ifelse(newx < 1, -1, ifelse(newx > side, -2, 0))
            withinboundsy = ifelse(newy < 1, -1, ifelse(newy > side, -2, 0))
            if withinboundsx < 0 || withinboundsy < 0
                return (withinboundsx, withinboundsy)
            else
                if notblocked
                    if rand(rng) < @inbounds shade_map[newx, newy] * prob_block
                        notblocked = false
                    end
                else
                    # if @inbounds (id = id_map[newx, newy]) != 0
                    if @inbounds farm_map[newx, newy] == 1
                        return (newx, newy)#, id
                    else
                        return (0,0)#, 0
                    end
                end
            end
        end
        advanced = false
        traveled += 0.5
    end

    if @inbounds farm_map[onx, ony] == 1
        (onx, ony)#, id
    else
        return (0,0)#, 0
    end
end

## Dispersal from outside the farm

function outside_spores!(model::SpatialRustABM)
    heading = model.current.wind_h
    expdist = Exponential(model.rustpars.wind_dst)
    outsp = model.outpour

    if isapprox(heading, 360.0; atol = 2.0) || isapprox(heading, 0.0; atol = 2.0)
        # cosd(2) ≈ 0.99939, which is just horizontal for a 100x100 farm
        tries = round(Int, outsp[1], RoundDown)
        deposited = fill((0,0), tries)
        i = 1
        while i <= tries
            deposited[i] = try_outside_disp!(model.rng, heading, model.farm_map, model.shade_map,
            model.rustpars, expdist, 1)
            i += 1
        end
    elseif heading < 90.0
        qs = [1,7,6]
        tries = round.(Int, outsp[qs], RoundDown)
        deposited = fill((0,0), sum(tries))
        i = 1
        ts = 0
        for (q, t) in zip(qs, tries)
            ts += t
            while i <= t
                deposited[i] = try_outside_disp!(model.rng, heading, model.farm_map, model.shade_map,
                model.rustpars, expdist, q)
                i += 1
            end
        end
    elseif isapprox(heading, 90.0; atol = 2.0) 
        tries = round(Int, outsp[6], RoundDown)
        deposited = fill((0,0), tries)
        i = 1
        while i <= tries
            deposited[i] = try_outside_disp!(model.rng, heading, model.farm_map, model.shade_map,
            model.rustpars, expdist, 6)
            i += 1
        end
    elseif heading < 180.0
        qs = [6,8,2]
        tries = round.(Int, outsp[qs], RoundDown)
        deposited = fill((0,0), sum(tries))
        i = 1
        ts = 0
        for (q, t) in zip(qs, tries)
            ts += t
            while i <= ts
                deposited[i] = try_outside_disp!(model.rng, heading, model.farm_map, model.shade_map,
                model.rustpars, expdist, q)
                i += 1
            end
        end
    elseif isapprox(heading, 180.0; atol = 2.0)
        tries = round(Int, outsp[2], RoundDown)
        deposited = fill((0,0), tries)
        i = 1
        while i <= tries
            deposited[i] = try_outside_disp!(model.rng, heading, model.farm_map, model.shade_map,
            model.rustpars, expdist, 2)
            i += 1
        end
    elseif heading < 270.0
        qs = [2,5,3]
        tries = round.(Int, outsp[qs], RoundDown)
        deposited = fill((0,0), sum(tries))
        i = 1
        ts = 0
        for (q, t) in zip(qs, tries)
            ts += t
            while i <= ts
                deposited[i] = try_outside_disp!(model.rng, heading, model.farm_map, model.shade_map,
                model.rustpars, expdist, q)
                i += 1
            end
        end
    elseif isapprox(heading, 270.0; atol = 2.0)
        tries = round(Int, outsp[3], RoundDown)
        deposited = fill((0,0), tries)
        i = 1
        while i <= tries
            deposited[i] = try_outside_disp!(model.rng, heading, model.farm_map, model.shade_map,
            model.rustpars, expdist, 3)
            i += 1
        end
    else
        qs = [3,4,1]
        tries = round.(Int, outsp[qs], RoundDown)
        deposited = fill((0,0), sum(tries))
        i = 1
        ts = 0
        for (q, t) in zip(qs, tries)
            ts += t
            while i <= ts
                deposited[i] = try_outside_disp!(model.rng, heading, model.farm_map, model.shade_map,
                model.rustpars, expdist, q)
                i += 1
            end
        end
    end

    for dep in filter!(t -> any(t .> 0), deposited)
        c = (model[id_in_position(dep, model)])
        if c.exh_countdown == 0
            c.newdeps += 1.0
            c.rusted = true
        end
    end
end

function try_outside_disp!(rng, heading::Float64, farm_map::Array{Int},
    shade_map::Array{Float64}, rp::RustPars, expdist::Exponential{Float64}, q::Int)

    startpos = starting_pos(rng, rp.map_side, q)
    distance = rand(rng, expdist) * (1 + rp.diff_wind)
    return gust(rng, startpos, distance, (heading + (rand(rng) * 30.0) - 15.0), farm_map, shade_map, rp)#[2]
end

function starting_pos(rng, side::Int, q::Int)
    if q == 1 # quadrant to the left
        return (rand(rng, 1:side), 1)
    elseif q == 7 # quadrant in the down-left diagonal
        quarter = fld(side, 4)
        randcoor = rand(rng, [1,2])
        if randcoor == 1
            return (rand(rng, (3*quarter+1):side), side)
        else
            return (1, rand(rng, 1:quarter))
        end
    elseif q == 6 # quadrant below
        return (side, rand(rng, 1:side))
    elseif q == 8 # quadrant in the down-right diagonal
        quarter = fld(side, 4)
        return (rand(rng, (3*quarter+1):side), side)[shuffle!(rng, [1,2])]
    elseif q == 2 # quadrant to the right
        return (rand(rng, 1:side), side)
    elseif q == 5 # quadrant in the up-right diagonal
        quarter = fld(side, 4)
        randcoor = rand(rng, [1,2])
        if randcoor == 1
            return (rand(rng, 1:quarter), side)
        else
            return (1, rand(rng, (3*quarter+1):side))
        end
    elseif q == 3 # quadrant above
        return (1, rand(rng, 1:side))
    else # q = 4, quadrant in the up-left diagonal
        quarter = fld(side, 4)
        return (rand(rng, 1:quarter), 1)[shuffle!(rng, [1,2])]
    end
end

function reintroduce_rusts!(model::SpatialRustABM, n_rusts::Int)
    activecofs = filter(active, model.agents)
    # n_rusts = min(n_rusts, length(activecofs))
    rusted_cofs = sample(model.rng, activecofs, n_rusts, replace = false)
    nl_distr = Binomial(model.rustpars.max_lesions - 1, 0.05)

    for rusted in rusted_cofs
        deposited = 0.0
        nl = n_lesions = 1 + rand(model.rng, nl_distr)
        ages = rusted.ages
        areas = rusted.areas
        spores = rusted.spores

        for _ in 1:nl
            area = rand(model.rng) * 0.2
            age = round(Int, area * 100.0)
            if 0.01 < area
                push!(ages, age)
                push!(areas, area)
                push!(spores, false)
            else
                deposited += 1.0
                n_lesions -= 1
            end
        end

        rusted.rusted = true
        rusted.deposited = deposited
        rusted.n_lesions = n_lesions
    end
end
