export dummyrun_spatialrust, simplerun

function dummyrun_spatialrust(steps::Int = 200, side::Int = 60, maxlesions::Int = 25; kwargs...)
    model = init_spatialrust(steps = steps, map_side = side, max_lesions = maxlesions; kwargs...)

    a_df, m_df = run!(model, dummystep, step_model!, steps;
        adata = [(:n_lesions, emedian, rusted), (tot_area, emedian, rusted), (:production, mean)],
        mdata = [incidence, n_rusts])
    rename!(a_df, [:step, :n_lesion_med, :tot_area_med, :production_mean])

    return a_df, m_df, model
end

function step_n!(model::SpatialRustABM, n::Int)
    s = 0
    while s < n
        step_model!(model)
        s += 1
    end
    return s
end

function step_while!(model::SpatialRustABM, s::Int, n::Int)
    while s < n && model.current.withinbounds
        step_model!(model)
        s += 1
    end
    return s
end

function simplerun(steps::Int = 365, retmodel::Bool = false; kwargs...)
    model = init_spatialrust(steps = steps; kwargs...)

    df = runsimple!(model, steps)

	if retmodel
    	return df, model
    else
    	return df
    end
end

function runsimple!(model::SpatialRustABM, steps::Int)
    meanshade = mean(model.shade_map)
    allcofs = model.agents
    ncofs = length(allcofs)
    sporepct = model.rustpars.spore_pct

    df = DataFrame(dayn = Int[],
        veg = Float64[], storage = Float64[], production = Float64[],
        indshade = Float64[], mapshade = Float64[],
        nl = Float64[], sumarea = Float64[], sumspore = Float64[],
        sporearea = Float64[],
        active = Float64[], farmprod = Float64[],
        nrusts = Int[], mages = Float64[], 
        occup = Float64[], incidence = Float64[],
    )
    for c in eachcol(df)
        sizehint!(c, steps)
    end

    s = 0
    while s < steps
        indshade = model.current.ind_shade

        sumareas = Iterators.filter(>(0.0), map(r -> sum(r.areas), allcofs))
        if isempty(sumareas)
            msuma = 0.0
            msumsp = 0.0
            mages = 0.0
        else
            msuma = mean(sumareas)
            sumspores = map(r -> sum(r.spores), allcofs)
            msumsp = mean(sumspores)
            mages = mean(map(a -> meanage(a.ages), allcofs))
        end

        push!(df, [
            model.current.days,
            mean(map(a -> a.veg, allcofs)),
            mean(map(a -> a.storage, allcofs)),
            mean(map(a -> a.production, allcofs)),
            indshade,
            indshade * meanshade,
            mean(map(a -> a.n_lesions, allcofs)),
            msuma,
            msumsp,
            mean(map(sporear, allcofs)) * sporepct,
            sum(map(active, allcofs)) / ncofs,
            copy(model.current.prod),
            sum(map(a -> a.rusted, allcofs)),
            mages,
            mean(a -> sum(a.n_lesions), allcofs),
            sum(a -> (a.n_lesions > 0), allcofs) / ncofs
        ])
        step!(model, dummystep, step_model!, 1)
        s += 1
    end

    indshade = model.current.ind_shade
    
    sumareas = Iterators.filter(>(0.0), map(r -> sum(r.areas), allcofs))
    if isempty(sumareas)
        msuma = 0.0
        msumsp = 0.0
        mages = 0.0
    else
        msuma = mean(sumareas)
        sumspores = map(r -> sum(r.spores), allcofs)
        msumsp = mean(sumspores)
        mages = mean(map(a -> meanage(a.ages), allcofs))
    end

    push!(df, [
        model.current.days,
        mean(map(a -> a.veg, allcofs)),
        mean(map(a -> a.storage, allcofs)),
        mean(map(a -> a.production, allcofs)),
        indshade,
        indshade * meanshade,
        mean(map(a -> a.n_lesions, allcofs)),
        msuma,
        msumsp,
        mean(map(sporear, allcofs)) * sporepct,
        sum(map(active, allcofs)) / ncofs,
        copy(model.current.prod),
        sum(map(a -> a.rusted, allcofs)),
        mages,
        mean(a -> sum(a.n_lesions) / 25, allcofs),
        sum(a -> (a.n_lesions > 0), allcofs) / ncofs
    ])

    return df
end

function sporear(a::Coffee)
    return sum((ar * sp for (ar,sp) in zip(a.areas, a.spores)), init = 0.0)
end

function meanage(ages)
    if isempty(ages)
        return 0.0
    else
        return mean(ages)
    end
end

emptymean(c) = isempty(c.areas) ? 0.0 : mean(filter(>(0.0), c.areas))

function cure(c::Coffee)
    c.rusted = false
    c.newdeps = 0.0
    c.deposited = 0.0
    c.n_lesions = 0
    empty!(c.areas)
    empty!(c.spores)
    empty!(c.ages)
end