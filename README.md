[![DOI](https://zenodo.org/badge/677233131.svg)](https://zenodo.org/badge/latestdoi/677233131)

# SpatialRust
 
A spatially explicit, process-based model of Coffee Leaf Rust (CLR) epidemics in coffee agroforestry systems. Built as the centerpiece of a doctoral dissertation in Biological Design (Arizona State University, 2023).
 
SpatialRust couples three interacting subsystems — shade tree growth and pruning dynamics, coffee plant physiology and biennial production cycles, and CLR infection, growth, sporulation, and dispersal — on a 100×100 spatial agent grid with a daily timestep. It is designed not to produce quantitative predictions, but to establish causal links between management decisions and epidemic and productivity outcomes over multi-year horizons.
 
## What the model does
 
- Simulates CLR epidemic dynamics under conventional and agroforestry-based farm management
- Encodes shade tree effects on local microclimate (temperature, radiation) and pathogen dispersal
- Models coffee plant photosynthesis, storage, and biennial production as a function of shading and disease load
- Allows exploration of management strategies: shade tree density, pruning schedules, spatial farm arrangement, fungicide timing, and inspection frequency
 
## Repository structure
 
```
src/
  SpatialRust.jl       # Package entry point
  ABM/                 # Core agent-based model
    CreateABM.jl       # Model initialization
    MainSetup.jl       # Parameter and environment setup
    MainStep.jl        # Daily simulation step
    FarmMap.jl         # Spatial farm layout
    ShadeMap.jl        # Shade tree spatial distribution
    CoffeeSteps.jl     # Coffee plant physiology and production
    CGrowerSteps.jl    # Coffee grower / management actions
    ShadeSteps.jl      # Shade tree growth and pruning dynamics
    RustGrowth.jl      # CLR lesion growth and sporulation
    RustDispersal.jl   # Wind and rain dispersal kernels
  QuickMetrics.jl      # Summary statistics and output metrics
  QuickRuns.jl         # Convenience functions for running scenarios
scripts/
  install.jl           # Dependency installation
  samplerun.jl         # Minimal working example
```
 
This repository contains the core simulation model. The full research pipeline — including the Approximate Bayesian Computation calibration (1 million parameter combinations parallelized across SLURM array jobs) and genetic algorithm optimization of management strategies — was developed on top of this model as part of the dissertation work.
 
## Usage
 
Install dependencies:
 
```
$ julia scripts/install.jl
```
 
Run a sample simulation:
 
```
$ julia scripts/samplerun.jl
```
 
Outputs are written to the `results/` folder. Parameter documentation is available in the [dissertation](https://hdl.handle.net/2286/R.2.N.190702).
 
## Citation
 
Vanegas Ferro, M. (2023). SpatialRust (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.8237935
 
## License
 
MIT. See LICENSE for details.
