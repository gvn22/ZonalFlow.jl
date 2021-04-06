module ZonalFlow

using OrdinaryDiffEq
using StochasticDiffEq
using DiffEqNoiseProcess
using DiffEqCallbacks
using RecursiveArrayTools
using FFTW
using LinearAlgebra
using Random,Distributions
using NPZ

include("structures.jl")
include("coefficients.jl")
include("equations.jl")
include("ics.jl")
include("noise.jl")
include("tools.jl")
include("solvers.jl")
include("analysis.jl")
include("writers.jl")

export Coefficients,
    Domain,
    PointJet,
    Kolmogorov,
    Stochastic,
    BetaPlane,
    NL,
    GQL,
    GCE2,
    CE2
    
# Solvers
export nl,gql,gce2,ce2

# ICs
export ic_pert_eqm,ic_eqm,ic_rand

# Analysis
export inversefourier
export energy,zonalenergy,fourierenergy,energyspectrum
export meanvorticity,zonalvelocity,meanzonalvelocity
export zonostrophy,energyinjectionrate
export adjacency

# Writers
export dumpenergy,dumpfields,dumpstats,dumpadjacency

end
