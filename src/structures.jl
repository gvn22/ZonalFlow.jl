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
    Δ::T
    τ::T
    isotropic::Bool
end
Stochastic(;kf,dk,ε,Δ=0.2,τ=0.0,isotropic=false) = Stochastic(kf,dk,ε,Δ,τ,isotropic)

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
eftype(::Type{<:AbstractProblem{T,F}}) where {T<:AbstractFloat,F<:AbstractForcing} = F

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

struct CE2 <: DSS
    poscheck::Bool
    poscheckat::Int
    eigmax::Bool
end
CE2(;poscheck=false,poscheckat=20,eigmax=false) = CE2(poscheck,poscheckat,eigmax)

struct GCE2 <: GSS
    Λ::Int
    poscheck::Bool
    poscheckat::Int
end
GCE2(Λ) = GCE2(Λ,false,1)
GCE2(;Λ,poscheck=false,poscheckat=20) = GCE2(Λ,poscheck,poscheckat)

label(eqs::NL) = "nl"
label(eqs::GQL) = "gql_$(eqs.Λ)"
label(eqs::CE2) = "ce2"
label(eqs::GCE2) = "gce2_$(eqs.Λ)"
label(eqs::Vector{AbstractEquations}) = [label(eq) for eq in eqs]

lambda(prob,eqs::NL)::Int = prob.d.nx-1
lambda(prob,eqs::GQL)::Int = eqs.Λ
lambda(prob,eqs::GCE2)::Int = eqs.Λ
lambda(prob,eqs::CE2)::Int = 0

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

"""
    Utilities for setting up ICs
"""
Base.similar(::Type{DNSField{T}},d::Domain{T};Λ::Int=d.nx-1) where T = similar(Field{T},2d.ny-1,Λ+1)

Base.similar(::Type{FirstCumulant{T}},d::Domain{T}) where T = similar(FirstCumulant{T},2d.ny-1)
Base.similar(::Type{SecondCumulant{T}},d::Domain{T}) where T = similar(SecondCumulant{T},2d.ny-1,2d.ny-1,d.nx-1)
Base.similar(::Type{DSSField{T}},d::Domain{T}) where T = ArrayPartition(similar(FirstCumulant{T},d),similar(SecondCumulant{T},d))

Base.similar(::Type{FieldBilinear{T}},d::Domain{T};Λ::Int) where T = similar(FieldBilinear{T},2d.ny-1,d.nx-Λ,2d.ny-1,d.nx-Λ)
Base.similar(::Type{GSSField{T}},d::Domain{T};Λ::Int) where T = ArrayPartition(similar(Field{T},d,Λ=Λ),similar(FieldBilinear{T},d,Λ=Λ))

Base.similar(eqs::Union{NL,GQL},d::Domain{T}) where T = similar(DNSField{T},d)
Base.similar(eqs::CE2,d::Domain{T}) where T = similar(DSSField{T},d)
Base.similar(eqs::GCE2,d::Domain{T}) where T = similar(GSSField{T},d,Λ=eqs.Λ)

Base.convert(::Union{NL,GQL},x::DNSField{T},d::Domain{T}) where T = x

function Base.convert(::CE2,x::DNSField{T},d::Domain{T}) where T
    (nx,ny) = size(d) #
    c1 = x[:,1]
    c2 = zeros(Complex{T},2ny-1,2ny-1,nx-1)
    @inbounds for m1=1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for n2=-ny+1:ny-1
                c2[n2+ny,n1+ny,m1] = conj(x[n2+ny,m1+1])*x[n1+ny,m1+1] # conj(2)*1
            end
        end
    end
    ArrayPartition(c1,c2)
end

function Base.convert(eqs::GCE2,x::DNSField{T},d::Domain{T}) where T
    (nx,ny),Λ = size(d),eqs.Λ
    c1 = x[:,1:Λ+1]
    c2 = zeros(Complex{T},2ny-1,nx-Λ,2ny-1,nx-Λ)
    @inbounds for m1=Λ+1:nx-1
        @inbounds for n1=-ny+1:ny-1
            @inbounds for m2=Λ+1:nx-1
                @inbounds for n2=-ny+1:ny-1
                    c2[n2+ny,m2-Λ,n1+ny,m1-Λ] = conj(x[n2+ny,m2+1])*x[n1+ny,m1+1] # conj(2)*1
                end
            end
        end
    end
    ArrayPartition(c1,c2)
end

function Base.convert(::CE2,x::GSSField{T},d::Domain{T}) where T
    (nx,ny),Λ = size(d),length(x.x[1][1,:])-1
    c1 = x.x[1][:,1]
    c2 = zeros(Complex{T},2ny-1,2ny-1,nx-1)
    @inbounds for n1=-ny+1:ny-1
        @inbounds for n2=-ny+1:ny-1
            @inbounds for m1=1:Λ
                c2[n2+ny,n1+ny,m1] = conj(x.x[1][n2+ny,m1+1])*x.x[1][n1+ny,m1+1]
            end
            @inbounds for m1=Λ+1:nx-1
                c2[n2+ny,n1+ny,m1] = x.x[2][n2+ny,m1-Λ,n1+ny,m1-Λ]
            end
        end
    end
    ArrayPartition(c1,c2)
end

abstract type AbstractParams{T <: AbstractFloat} end

struct NLParams{T} <: AbstractParams{T}
    nx::Int
    ny::Int
    A::Array{Complex{T},1}
    B::Array{Complex{T},2}
    C⁺::Array{T,4}
    C⁻::Array{T,4}
    F::Array{T,2}
end

struct GQLParams{T} <: AbstractParams{T}
    nx::Int
    ny::Int
    Λ::Int
    A::Array{Complex{T},1}
    B::Array{Complex{T},2}
    C⁺::Array{T,4}
    C⁻::Array{T,4}
    F::Array{T,2}
end

struct CE2Params{T} <: AbstractParams{T}
    nx::Int
    ny::Int
    A::Array{Complex{T},1}
    B::Array{Complex{T},2}
    C⁺::Array{T,4}
    C⁻::Array{T,4}
    F::ArrayPartition{T,Tuple{Array{T,1},Array{T,3}}}
    dy::SecondCumulant{T}
    temp::SecondCumulant{T}
end
CE2Params(nx,ny,A::Array{Complex{T},1},B::Array{Complex{T},2},C⁺::Array{T,4},C⁻::Array{T,4},F::ArrayPartition{T,Tuple{Array{T,1},Array{T,3}}}) where T = CE2Params(nx,ny,A,B,C⁺,C⁻,F,zeros(Complex{T},2ny-1,2ny-1,nx-1),zeros(Complex{T},2ny-1,2ny-1,nx-1))

struct GCE2Params{T} <: AbstractParams{T}
    nx::Int
    ny::Int
    Λ::Int
    A::Array{Complex{T},1}
    B::Array{Complex{T},2}
    C⁺::Array{T,4}
    C⁻::Array{T,4}
    F::ArrayPartition{T,Tuple{Array{T,2},Array{T,4}}}
    dx::Array{Complex{T},2}
    dy::FieldBilinear{T}
    temp::FieldBilinear{T}
end
GCE2Params(nx,ny,Λ,A::Array{Complex{T},1},B::Array{Complex{T},2},C⁺::Array{T,4},C⁻::Array{T,4},F::ArrayPartition{T,Tuple{Array{T,2},Array{T,4}}}) where T = GCE2Params(nx,ny,Λ,A,B,C⁺,C⁻,F,zeros(Complex{T},2ny-1,Λ+1),zeros(Complex{T},2ny-1,nx-Λ,2ny-1,nx-Λ),zeros(Complex{T},2ny-1,nx-Λ,2ny-1,nx-Λ))

params(d,eqs::NL,p) = NLParams(d.nx,d.ny,p...)
params(d,eqs::GQL,p) = GQLParams(d.nx,d.ny,eqs.Λ,p...)
params(d,eqs::CE2,p) = CE2Params(d.nx,d.ny,p[1],p[2],p[3],p[4],p[5])
params(d,eqs::GCE2,p) = GCE2Params(d.nx,d.ny,eqs.Λ,p[1],p[2],p[3],p[4],p[5])
