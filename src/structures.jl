"""
Coefficient type
"""
abstract type AbstractCoefficients{T <: AbstractFloat} end

struct Coefficients{T} <: AbstractCoefficients{T}
    Ω::T
    θ::T
    μ::T
    ν::T
    ν₄::T
    linear::Bool
end
Coefficients(;Ω,θ,μ,ν,ν₄,linear=false) = Coefficients(Ω,θ,μ,ν,ν₄,linear)
Base.eltype(::Type{<:AbstractCoefficients{T}}) where {T <: AbstractFloat} = T

"""
Forcing types
"""
abstract type AbstractForcing{T <: AbstractFloat} end

struct PointJet{T} <: AbstractForcing{T}
    Ξ::T
    Δθ::T
    τ::T
end
PointJet(;Ξ,Δθ,τ) = PointJet(Ξ,Δθ,τ)

struct Stochastic{T}  <: AbstractForcing{T}
    kf::Int
    dk::Int
    ε::T
end
Stochastic(;kf,dk,ε) = Stochastic(kf,dk,ε)

struct Kolmogorov{T} <: AbstractForcing{T}
    A₁::T
    A₄::T
end
Kolmogorov(;A₁,A₄) = Kolmogorov(A₁,A₄)

"""
    Domain and problem types
"""
abstract type AbstractDomain{T <: AbstractFloat} end

struct Domain{T} <: AbstractDomain{T}
    lx::T
    ly::T
    nx::Int
    ny::Int
end
Domain(;extent,res) = Domain(extent[1],extent[2],res[1],res[2])
Base.length(d::Domain) = ((d.lx,d.ly))
Base.size(d::Domain) = ((d.nx,d.ny))

abstract type AbstractProblem{T<:AbstractFloat,F<:AbstractForcing} end

struct BetaPlane{T,F} <: AbstractProblem{T,F}
    d::Domain{T}
    c::Coefficients{T}
    f::F
end
BetaPlane(coeffs,forcing;extent,res) = BetaPlane(Domain(extent=extent,res=res),coeffs,forcing)

Base.eltype(::Type{<:AbstractProblem{T,F}}) where {T<:AbstractFloat,F<:AbstractForcing} = T
eftype(::Type{<:AbstractProblem{T,F}}) where {T,F} = F

"""
    Fields and equations types
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
GCE2(Λ) = GCE2(Λ,false,1)
GCE2(;Λ,poscheckat=20,poscheck=false) = GCE2(Λ,poscheck,poscheckat)

const Field{T} = Array{Complex{T},2} where {T <: AbstractFloat}
const FirstCumulant{T} = Array{Complex{T},1} where {T <: AbstractFloat}
const SecondCumulant{T} = Array{Complex{T},3} where {T <: AbstractFloat}
const FieldBilinear{T} = Array{Complex{T},4} where {T <: AbstractFloat}

const DNSField{T} = Field{T} where {T <: AbstractFloat}
const DSSField{T} = ArrayPartition{Complex{T},Tuple{FirstCumulant{T},SecondCumulant{T}}} where {T <: AbstractFloat}
const GSSField{T} = ArrayPartition{Complex{T},Tuple{Field{T},FieldBilinear{T}}} where {T <: AbstractFloat}

function Base.getproperty(f::Union{DSSField{T},GSSField{T}},v::Symbol) where {T<:AbstractFloat}
    v == :l && return getfield(f,:x)[1]
    v == :h && return getfield(f,:x)[2]
    return getfield(f,v)
end
