module SpatialRust

using Agents, DataFrames, Distributions, Random
using StatsBase: sample, weights

include("ABM/MainSetup.jl")

const SpatialRustABM = Agents.SingleContainerABM{
    Agents.GridSpaceSingle{2, false}, Coffee, Vector{Coffee},
    typeof(Agents.Schedulers.fastest), SpatialRust.Props, Random.Xoshiro
}

include("ABM/CreateABM.jl")
include("ABM/FarmMap.jl")
include("ABM/ShadeMap.jl")

include("ABM/MainStep.jl")
include("ABM/ShadeSteps.jl")
include("ABM/CoffeeSteps.jl")
include("ABM/RustGrowth.jl")
include("ABM/RustDispersal.jl")
include("ABM/CGrowerSteps.jl")

include("QuickRuns.jl")
include("QuickMetrics.jl")

export SpatialRustABM

end
