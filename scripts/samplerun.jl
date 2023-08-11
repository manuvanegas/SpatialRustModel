# Load Packages
import Pkg
Pkg.activate(".") 
using CSV, DataFrames, SpatialRust

# Run one simulation using default parameters
df = simplerun(365, false)

# If you want to get the model object as well, run
# df, model = simplerun(365, true)

# Modify 365 to the number of days you want to simulate
# The rest of the model's input parameters can be modified via keywords.
# Example: to run the model using a custom shade pruning schedule, run
# df = simplerun(365, prune_sch = [15, 196, 285])
# Or to run it with this custom schedule and a given mean temperature
# df = simplerun(365, prune_sch = [15, 196, 285], mean_temp = 19.0)
# A complete list of input parameters and their default values will be made available soon

# Write results
mkpath("results/")
CSV.write("results/samplerun.csv", df)