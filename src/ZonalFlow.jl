module ZonalFlow

using OrdinaryDiffEq
using StochasticDiffEq
using DiffEqCallbacks
using RecursiveArrayTools
using FFTW
using LinearAlgebra
using Random
using Distributions
using NPZ
using JLD2

include("structures.jl")
include("ic.jl")
include("coeffs.jl")
include("solve.jl")
include("equations.jl")
include("tools.jl")
include("analysis.jl")
include("writers.jl")
# include("deprecated.jl")
# include("future.jl")

export  Domain,
        Coefficients,
        PointJet,
        Kolmogorov,
        Stochastic,
        BetaPlane

export  NL,
        GQL,
        GCE2,
        CE2

export  integrate

end
