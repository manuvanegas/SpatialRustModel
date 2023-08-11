export Coffee, init_spatialrust, createweather, create_farm_map, step_n!

# Coffee agent type
@agent Coffee GridAgent{2} begin
    sunlight::Float64 # let through by shade trees
    veg::Float64
    storage::Float64
    production::Float64
    exh_countdown::Int
    rust_gr::Float64

    #Rust
    rusted::Bool
    newdeps::Float64
    deposited::Float64 
    n_lesions::Int
    ages::Vector{Int}
    areas::Vector{Float64}
    spores::Vector{Bool}
end

# Coffee constructor
function coffee(id, pos, max_lesions::Int, rust_gr::Float64;
    sunlight::Float64 = 1.0, veg::Float64 = 2.0, storage::Float64 = 100.0)

    Coffee(
        id, pos, sunlight, veg, storage, 0.0, 0, rust_gr,
        false, 0.0, 0.0, 0,
        sizehint!(Int[], max_lesions), sizehint!(Float64[], max_lesions), sizehint!(Bool[], max_lesions),
    )
end

function init_spatialrust(;
    seed::Int = 0,
    start_days_at::Int = 0,
    steps::Int = 500,                       # simulation steps
    p_rusts::Float64 = 0.01,                # % of initial rusts (# of initial clusters, if > 1)
    p_row::Int = 0,                         # parameter combination number (for ABC)
    rep::Int = 0,                           # repetition number (for other exps)

    # weather parameters
    rain_prob::Float64 = 0.8,
    wind_prob::Float64 = 0.7,
    mean_temp::Float64 = 22.0,
    rain_data::Vector{Bool} = Bool[],       # if provided, rain_prob is ignored
    wind_data::Vector{Bool} = Bool[],       # if provided, wind_prob is ignored
    temp_data::Vector{Float64} = Float64[], # if provided, mean_temp is ignored

    # coffee parameters
    veg_d::Int = 1,                         # start of vegetative growth   
    rep_d::Int = 135,                       # start of reproductive growth
    f_avail::Float64 = 0.5,                 # fraction of daily assimilates available for allocation
    phs_max::Float64 = 0.2,                 # maximum assimilation rate
    k_sl::Float64 = 0.05,                   # half-rate constant for sunlight-dependent photosynthesic rate  
    k_v::Float64 = 0.2,                     # half-rate constant for veg-dependent photosynthesic rate 
    photo_frac::Float64 = 0.2,              # fraction of veg tissue that is photosynthetic
    phs_veg::Float64 = 0.8,                 # fraction of photosynthetate allocated and converted to veg 
    μ_veg::Float64 = 0.01,                  # rate of veg biomass loss
    phs_sto::Float64 = 0.6,                 # fraction of photosynthetate allocated and converted to storage
    res_commit::Float64 = 0.25,             # scaling par to determine resources commited to production
    μ_prod::Float64 = 0.01,                 # rate of production biomass loss

    # rust parameters
    max_lesions::Int64 = 25,                # maximum number of rust lesions
    rain_washoff::Float64 = 0.25,           # rain wash-off " " " "; Avelino et al., 2020
    viab_loss::Float64 = 0.75,              # Deposited spore viability loss
    fung_inf::Float64 = 0.9,                # infection prob under fungicide mod
    fung_gro_prev::Float64 = 0.3,           # fungicide mod to growth rate on preventive fungicide
    fung_gro_cur::Float64 = 0.6,            # fungicide mod to growth rate on curative fungicide
    fung_spor_prev::Float64 = 0.0,          # fungicide mod to spor prob on preventive fungicide
    fung_spor_cur::Float64 = 0.5,           # fungicide mod to spor prob on curative fungicide
    fung_reduce::Float64 = 0.95,            # fungicide's area reduction rate
    
    exh_countdown::Int = 731,               # days to count after plant has been exhausted (2-3 y to resume production) 

    map_side::Int = 100,                    # side size
    diff_splash::Float64 = 2.5,             # times rain distance due to enhanced kinetic e (shade) (Avelino et al. 2020: "Kinetic energy was twice as high"+ Gagliardi et al. 2021: TKE ~ 5xOpenness%)
    diff_wind::Float64 = 3.0,               # additional wind distance due to increased openness (Pezzopane

    # farm management
    harvest_day::Int = 365,                 # 365 or 182, depending on the region
    prune_sch::Vector{Int} = [74, 196,-1],
    inspect_period::Int = 7,                # days
    fungicide_sch::Vector{Int} = [91,273,-1],
    fung_stratg::Symbol = :cal,             # fungicide use strategy: :cal, :cal_incd, :incd, :flor
    incidence_thresh::Float64 = 0.1,        # incidence that triggers fungicide use (10% from Cenicafe's Boletin 36)
    max_fung_sprayings::Int = 3,            # maximum fung treatments per year

    prune_cost::Float64 = 1.656,            # per shade
    inspect_cost::Float64 = 0.006,          # per coffee inspected
    fung_cost::Float64 = 0.179,             # per coffee
    other_costs::Float64 = 0.088,           # other costs considered (weeding)
    fixed_costs::Float64 = 7600.6,          # production-independent costs
    coffee_price::Float64 = 1.0,

    post_prune::Vector{Float64} = [0.3, 0.5, 0.0], # individual shade level after pruning
    inspect_effort::Float64 = 0.01,         # % coffees inspected each time
    rm_lesions::Int = 2,                    # if rust is found, remove n lesions
    fung_effect::Int = 30,                  # length of fungicide effect

    # shade parameters
    max_shade::Float64 = 0.8,               # maximum individual shade
    shade_g_rate::Float64 = 0.015,           # shade growth rate
    shade_r::Int = 3,                       # radius of influence of shades

    # farm map
    farm_map::Array{Int} = Int[],           # if provided, parameters below are ignored
    common_map::Symbol = :none,             # :fullsun or :regshaded
    row_d::Int = 2,                         # distance between rows (options: 1, 2, 3)
    plant_d::Int = 1,                       # distance between plants (options: 1, 2)
    shade_d::Int = 6,                       # distance between shades (only considered when :regular)
    shade_pattern::Symbol = :regular,       # or :rand
    barrier_rows::Int = 2,                  # or 2 = double
    barriers::NTuple{2, Int} = (1, 0),      # barrier arrangement: 1->internal(0, 1, or 2),2->edges(0 or 1)

    # keeping calibrated pars here for quick reference
    # temp_cooling::Float64 = 1.0665,       # temp reduction due to shade
    # max_inf::Float64 = 1.0301,            # Max infection probability
    # light_inh::Float64 = 0.5168,          # UV inactivation prob under 100% sunlight 
    # rep_inf::Float64 = 0.2113,            # Weight of reprod growth on infection prob
    # rust_gr::Float64 = 0.0867,            # basic rust area growth rate
    # opt_temp::Float64 = 22.8716,          # optimal rust growth temp
    # rep_spo::Float64 = 0.3926,            # Effect on reprod growth on sporul
    # pdry_spo::Float64 = 0.7078,           # Prob of sporulation without rain
    # temp_ampl::Float64 = 5.0367,          # (max temp - optimal temp)
    # rep_gro::Float64 = 0.01,            # resource sink effect on area growth
    # spore_pct::Float64 = 0.3836,            # % of area that sporulates
    # rust_paras::Float64 = 0.0155,         # resources taken per unit of total area
    # rain_dst::Float64 = 2.2071,           # mean distance of spores dispersed by rain
    # tree_block::Float64 = 0.4889,         # prob a tree will block rust dispersal
    # wind_dst::Float64 = 5.7806,           # mean distance of spores dispersed by wind
    # shade_block::Float64 = 0.7898,        # prob of Shades blocking a wind dispersal event
    # les_surv::Float64 = 0.4338,           # proportion of lesions surviving to next cycle
)

    rng = seed == 0 ? Xoshiro() : Xoshiro(seed)

    w = noisedweather(rain_data, wind_data, temp_data, rain_prob, wind_prob, mean_temp, steps, rng)

    if isempty(farm_map)
        if common_map == :none
            farm_map = create_farm_map(map_side, row_d, plant_d, shade_d, shade_pattern, barrier_rows, barriers)
        elseif common_map == :fullsun
            farm_map = create_fullsun_farm_map(map_side, row_d, plant_d)
        elseif common_map == :regshaded
            farm_map = create_regshaded_farm_map(map_side, row_d, plant_d, shade_d)
        end
    else
        map_side = size(farm_map)[1]
    end

    smap = create_shade_map(farm_map, shade_r, map_side, common_map)

    cp = CoffeePars(
        veg_d, rep_d, f_avail * phs_max, k_sl, k_v, photo_frac,
        phs_veg, μ_veg, phs_sto, res_commit, μ_prod, exh_countdown
    )

    rp = RustPars(
        # Infection
        max_lesions, 1.0665, 0.5168, rain_washoff, 0.2113, viab_loss, 1.0301, 
        # Sporulation/Growth
        0.0867, 22.8716, 0.3926, 0.7078, -(1.0/5.0367^2), 0.01, 0.3836, 
        fung_inf, fung_gro_prev, fung_gro_cur, fung_spor_prev, fung_spor_cur, fung_reduce,
        # Parasitism
        0.0155, exh_countdown,
        # Dispersal
        map_side, 2.2071, diff_splash, 0.4889, 5.7806, diff_wind, 0.7898
    )

    pruneskept = filter!(i -> prune_sch[i] > 0, sortperm(prune_sch))
    prune_sch = prune_sch[pruneskept]
    post_prune = post_prune[pruneskept]
    if (l = length(prune_sch) == 2) && !allunique(prune_sch)
        prune_sch = keepat!(prune_sch, 1)
        post_prune = minimum(post_prune)
    elseif l == 3 && !allunique(prune_sch)
        if prune_sch[1] == prune_sch[3]
            prune_sch = keepat!(prune_sch, 1)
            post_prune = minimum(post_prune)
        elseif prune_sch[2] == prune_sch[3]
            pop!(prune_sch)
            ts = pop!(post_prune)
            post_prune[2] = min(post_prune[2], ts)
        elseif prune_sch[1] == prune_sch[2]
            popat!(prune_sch, 2)
            ts = popat!(post_prune, 2)
            post_prune[1] = min(post_prune[1], ts)
        end
    end
    prune_sch = Tuple(prune_sch)
    post_prune = Tuple(post_prune)
    if fung_stratg == :cal || fung_stratg == :cal_incd
        fungicide_sch = Tuple(sort!(filter!(>(0), fungicide_sch)))
    elseif fung_stratg == :flor
        fungicide_sch = (rep_d + 60, rep_d + 105)
    else
        fungicide_sch = ()
    end
    
    n_shades = count(farm_map .== 2)
    n_coffees = count(farm_map .== 1)

    n_inspect = trunc(Int, inspect_effort * n_coffees)
    inspcst = inspect_cost * (1.0 + 0.2 * (rm_lesions - 1))

    mp = MngPars{length(prune_sch),length(fungicide_sch)}(
        harvest_day, prune_sch,
        inspect_period, rm_lesions,
        fungicide_sch, fung_stratg, incidence_thresh, max_fung_sprayings,
        #
        n_shades, rand(rng, Normal(prune_cost, prune_cost * 0.05)) * n_shades,
        n_coffees, rand(rng, Normal(inspcst, inspcst * 0.05)) * n_inspect, rand(rng, Normal(fung_cost, fung_cost * 0.05)) * n_coffees,
        rand(rng, Normal(other_costs, other_costs * 0.05)), rand(rng, Normal(fixed_costs, fixed_costs * 0.05)), coffee_price,
        #
        0.4338, post_prune, n_inspect, fung_effect,
        max_shade, shade_g_rate, shade_r
    )

    doy = start_days_at == 0 ? veg_d - 1 : start_days_at

    b = Books(
        doy, 0, ind_shade_i(shade_g_rate, max_shade, doy, post_prune, prune_sch),
        0.0, false, false, 0.0, 0, 0, 0.0, 0.0, 0.0, true, true, 0.0
    )

    return init_abm_obj(Props(w, cp, rp, mp, b, farm_map, smap, zeros(8)), rng, p_rusts)
end

# Definitions of the different parameter structs

struct Weather
    rain_data::Vector{Bool}
    wind_data::Vector{Bool}
    temp_data::Vector{Float64}
end

struct CoffeePars
    veg_d::Int
    rep_d::Int
    photo_const::Float64 # f_avail * phs_max * photo_frac
    k_sl::Float64
    k_v::Float64
    photo_frac::Float64
    phs_veg::Float64
    μ_veg::Float64
    phs_sto::Float64
    res_commit::Float64
    μ_prod::Float64
    exh_countdown::Int          
end

struct RustPars
    # growth
    max_lesions::Int64
    temp_cooling::Float64
    light_inh::Float64
    rain_washoff::Float64
    rep_inf::Float64
    viab_loss::Float64
    max_inf::Float64
    rust_gr::Float64
    opt_temp::Float64
    rep_spo::Float64
    pdry_spo::Float64
    temp_ampl_c::Float64
    rep_gro::Float64
    spore_pct::Float64
    fung_inf::Float64 
    fung_gro_prev::Float64
    fung_gro_cur::Float64 
    fung_spor_prev::Float64
    fung_spor_cur::Float64
    fung_reduce::Float64
    # parasitism
    rust_paras::Float64
    exh_countdown::Int
    # dispersal
    map_side::Int
    rain_dst::Float64
    diff_splash::Float64
    tree_block::Float64
    wind_dst::Float64
    diff_wind::Float64
    shade_block::Float64
end

struct MngPars{N,M}
    # action scheduling
    harvest_day::Int
    prune_sch::NTuple{N,Int}
    inspect_period::Int
    rm_lesions::Int
    fungicide_sch::NTuple{M,Int}
    fung_stratg::Symbol
    incidence_thresh::Float64
    max_fung_sprayings::Int
    # financials
    n_shades::Int
    tot_prune_cost::Float64
    n_cofs::Int
    tot_inspect_cost::Float64
    tot_fung_cost::Float64
    other_costs::Float64
    fixed_costs::Float64
    coffee_price::Float64
    # others
    lesion_survive::Float64
    # les_surv::Float64
    post_prune::NTuple{N, Float64}
    n_inspected::Int
    fung_effect::Int
    # by_fragments::Bool = true,            # apply fungicide differentially by fragments? (possible future development)
    max_shade::Float64
    shade_g_rate::Float64
    shade_r::Int
end

mutable struct Books
    days::Int                               # same as ticks unless start_days_at != 0
    ticks::Int                              # initialized as 0 but the 1st thing that happens is +=1, so it effectvly starts at 1
    ind_shade::Float64
    temperature:: Float64
    rain::Bool
    wind::Bool
    wind_h::Float64
    fungicide::Int
    fung_count::Int
    obs_incidence::Float64
    costs::Float64
    prod::Float64
    inbusiness::Bool
    withinbounds::Bool
    shadeacc::Float64
end

struct Props
    weather::Weather
    coffeepars::CoffeePars
    rustpars::RustPars
    mngpars::MngPars
    current::Books
    farm_map::Array{Int, 2}
    shade_map::Array{Float64}
    outpour::Vector{Float64}
    # 8 positions, one per direction the spores can leave the farm (from (0,0) which is the farm)
    # indexing is weird (also remember, julia goes column-first): 
    # 1 -> (0,-1), 2 -> (0,1), 3 -> (-1,0), 4 -> (-1,-1), 5 -> (-1,1), 6 -> (1,0), 7 -> (1,-1), 8 ->(1,1)
end

# Other fncs

function createweather(rain_prob, wind_prob, mean_temp, steps, rng)
    raindata = rand(rng, steps) .< rain_prob
    return Weather(
        raindata,
        rand(rng, steps) .< ifelse.(raindata, wind_prob - 0.1, wind_prob + 0.1),
        round.(rand(rng, Normal(mean_temp, 0.75), steps), digits = 2)
    )
end

function noisedweather(rain_data, wind_data, temp_data, rain_prob, wind_prob, mean_temp, steps, rng)
    raindata = isempty(rain_data) ? rand(rng, steps) .< rain_prob : addnoise(rain_data, steps, rng)
    wpdn = wind_prob - 0.1
    wpup = wind_prob + 0.1
    return Weather(
        raindata,
        isempty(wind_data) ? rand(rng, steps) .< ifelse.(raindata, wpdn, wpup) : addnoise(wind_data, steps, rng),
        isempty(temp_data) ? rand(rng, Normal(mean_temp, 0.5), steps) : addqnoise(temp_data, steps, rng)
    )
end

function addnoise(v::Vector{Bool}, steps, rng)
    @inbounds for (i, b) in enumerate(v)
        rand(rng) < 0.05 && (v[i] = !b)
    end
    return resize!(v, steps)
end

function addqnoise(v::Vector{Float64}, steps, rng)
    resize!(v, steps)
    return v .+ rand(rng, Normal(0, 0.05), steps)
end

# Calculate initial ind_shade
function ind_shade_i(
    shade_g_rate::Float64,
    max_shade::Float64,
    start_day_at::Int,
    post_prune::NTuple{N, Float64},
    prune_sch::NTuple{N, Int}) where {N}

    if isempty(prune_sch)
        return max_shade
    else
        day = start_day_at > 0 ? start_day_at : 1
        # calculate elapsed days since last prune
        prune_diff = filter(>(0), day .- prune_sch)
        if isempty(prune_diff)
            lastprune, prune_i = findmax(prune_sch)
            last_prune = 365 + day - lastprune
            pruned_to = post_prune[prune_i]
        else
            lastprune, prune_i = findmin(prune_diff)
            last_prune = minimum(prune_diff)
            pruned_to = post_prune[prune_i]
        end
        # logistic equation to determine starting shade level
        return max_shade * (pruned_to / (pruned_to + (max_shade - pruned_to) * exp(-(shade_g_rate * last_prune))))
    end
end
