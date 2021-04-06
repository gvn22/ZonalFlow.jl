"""
Coefficient type
"""

abstract type AbstractCoefficients end

struct Coefficients{T} <: AbstractCoefficients where {T <: AbstractFloat}
    Ω::T
    θ::T
    μ::T
    ν::T
    ν₄::T
    linear::Bool
end
Coefficients(;Ω,θ,μ,ν,ν₄,linear=false) = Coefficients(Ω,θ,μ,ν,ν₄,linear)

"""
Forcing types
"""

abstract type AbstractForcing end

struct PointJet{T} <:AbstractForcing where {T <: AbstractFloat}
    Ξ::T
    Δθ::T
    τ::T
end
PointJet(;Ξ,Δθ,τ) = PointJet(Ξ,Δθ,τ)

struct Stochastic{T}  <: AbstractForcing where {T <: AbstractFloat}
    kf::Int
    dk::Int
    ε::T
end
Stochastic(;kf,dk,ε) = Stochastic(kf,dk,ε)

struct Kolmogorov{T} <: AbstractForcing where {T <: AbstractFloat}
    A₁::T
    A₄::T
end
Kolmogorov(;A₁,A₄) = Kolmogorov(A₁,A₄)

"""
Domain and Model types
"""

abstract type AbstractDomain end

struct Domain{T} <: AbstractDomain where {T <: AbstractFloat}
    lx::T
    ly::T
    nx::Int
    ny::Int
end
Domain(;extent,res) = Domain(extent[1],extent[2],res[1],res[2])
Base.length(d::Domain) = ((d.lx,d.ly))
Base.size(d::Domain) = ((d.nx,d.ny))

abstract type AbstractProblem end

struct BetaPlane <: AbstractProblem
    d::Domain
    c::Coefficients
    f::AbstractForcing
end
BetaPlane(dom::Domain,coeffs::AbstractCoefficients,forcing::AbstractForcing) = BetaPlane(dom,coeffs,forcing)
BetaPlane(coeffs::AbstractCoefficients,forcing::AbstractForcing;extent,res) = BetaPlane(Domain(extent=extent,res=res),coeffs,forcing)

"""
Solver types
"""

abstract type AbstractEquations end

abstract type DNS <: AbstractEquations end
abstract type DSS <: AbstractEquations end
abstract type GSS <: AbstractEquations end

struct NL <: DNS end

struct GQL <: DNS
    Λ::Int
end
GQL(;Λ) = GQL(Λ)

struct CE2 <: DSS end

struct GCE2 <: GSS
    Λ::Int
    poscheck::Bool
    poscheckat::Int
end
GCE2(;Λ,poscheckat=20,poscheck=false) = GCE2(Λ,poscheck,poscheckat)

"""
"""
const Field{T} = Array{Complex{T},2} where {T <: AbstractFloat}
const FirstCumulant{T} = Array{Complex{T},1} where {T <: AbstractFloat}
const SecondCumulant{T} = Array{Complex{T},3} where {T <: AbstractFloat}
const FieldBilinear{T} = Array{Complex{T},4} where {T <: AbstractFloat}

const DNSField{T} = Field{T} where {T <: AbstractFloat}
const DSSField{T} = ArrayPartition{Complex{T},Tuple{FirstCumulant{T},SecondCumulant{T}}} where {T <: AbstractFloat}
const GSSField{T} = ArrayPartition{Complex{T},Tuple{Field{T},FieldBilinear{T}}} where {T <: AbstractFloat}
