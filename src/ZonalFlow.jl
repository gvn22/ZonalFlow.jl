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

include("coefficients.jl")
include("equations.jl")
include("ics.jl")
include("noise.jl")
include("tools.jl")
include("solvers.jl")
include("analysis.jl")
include("writers.jl")

export nl,gql,gce2,ce2 # solvers
export ic_pert_eqm,ic_eqm,ic_rand # ics
export energy,zonalenergy,inversefourier,meanvorticity,zonalvelocity,meanzonalvelocity,fourierenergy,zonostrophy,injectionrate # analysis
export dumpenergy,dumpfields,dumpstats

end
