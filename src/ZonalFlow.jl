module ZonalFlow

using OrdinaryDiffEq
using DiffEqCallbacks
using RecursiveArrayTools
using FFTW
using LinearAlgebra

include("coefficients.jl")
include("equations.jl")
include("ics.jl")
include("tools.jl")
include("solvers.jl")
include("analysis.jl")

export nl,gql,gce2,ic_pert_eqm,ic_eqm
export energy,zonalenergy,inversefourier,meanvorticity,zonalvelocity,fourierenergy

end
