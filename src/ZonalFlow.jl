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

"""
    Field type aliases
"""
const Field{T} = Array{Complex{T},2} where T
const FirstCumulant{T} = Array{Complex{T},1} where T
const SecondCumulant{T} = Array{Complex{T},3} where T
const FieldBilinear{T} = Array{Complex{T},4} where T

const DNSField{T} = Field{T} where T
const DSSField{T} = ArrayPartition{Complex{T},Tuple{FirstCumulant{T},SecondCumulant{T}}} where T
const GSSField{T} = ArrayPartition{Complex{T},Tuple{Field{T},FieldBilinear{T}}} where T

include("coefficients.jl")
include("equations.jl")
include("ics.jl")
include("noise.jl")
include("tools.jl")
include("solvers.jl")
include("analysis.jl")
include("writers.jl")

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
