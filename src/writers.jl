"""
    Overloads for Base.write
"""
Base.write(prob,eqs::Vector{AbstractEquations},sols;dn::String,labels::Vector{String}=label(eqs)) = foreach(x->write(prob,x[1],x[2],dn=dn,fn=x[3]),zip(eqs,sols,labels))

function Base.write(prob,eqs,sol;dn::String="",fn::String)

    (lx,ly),(nx,ny),Λ = length(prob.d),size(prob.d),lambda(prob,eqs)
    t,u = sol.t,sol.u

    E = dumpenergy(lx,ly,nx,ny,t,u)
    F = typeof(eqs) == GCE2 ? dumpfields(lx,ly,nx,ny,t,u,Λ=Λ) : dumpfields(lx,ly,nx,ny,t,u)
    S = typeof(eqs) == CE2 ? dumpstats(prob,u) : Dict("empty"=>0)

    mkpath(dn)
    NPZ.npzwrite(fn*".npz",merge(Dict("t"=>sol.t),E,F,S))
end

function dumpenergy(lx::T,ly::T,nx::Int,ny::Int,t::Array{T,1},u;t0::Float64=200.0) where {T <: AbstractFloat}
    Et,Zt = energy(lx,ly,nx,ny,u)
    Etav,Ztav = energy(lx,ly,nx,ny,t,u,t0=t0)
    Emt,Zmt = zonalenergy(lx,ly,nx,ny,u)
    Emtav,Zmtav = zonalenergy(lx,ly,nx,ny,t,u,t0=t0)
    Dict("Zt"=>Zt,"Ztav"=>Ztav,"Et"=>Et,"Etav"=>Etav,"Emt"=>Emt,"Emtav"=>Emtav)
end

function dumpfields(lx::T,ly::T,nx::Int,ny::Int,t::Array{Float64,1},u::Array{DNSField{T},1};Λ::Int=nx-1) where {T <: AbstractFloat}
    Emn = energyspectrum(lx,ly,nx,ny,u,Λ=Λ)
    Vxy = vorticity(nx,ny,u,Λ=Λ)
    Uxy = zonalvelocity(lx,ly,nx,ny,u,Λ=Λ)
    Vyt = zonalvorticity(nx,ny,u)
    Dict("Emn"=>Emn,"Vxy"=>Vxy,"Uxy"=>Uxy,"Vyt"=>Vyt)
end

function dumpfields(lx::T,ly::T,nx::Int,ny::Int,t::Array{Float64,1},u::Array{DSSField{T},1}) where {T <: AbstractFloat}
    Emn = energyspectrum(lx,ly,nx,ny,u)
    Vxy = vorticity(nx,ny,u)
    Uxy = zonalvelocity(lx,ly,nx,ny,u)
    Vyt = zonalvorticity(nx,ny,u)
    Dict("Emn"=>Emn,"Vxy"=>Vxy,"Uxy"=>Uxy,"Vyt"=>Vyt)
end

function dumpfields(lx::T,ly::T,nx::Int,ny::Int,t::Array{Float64,1},u::Array{GSSField{T},1};Λ::Int) where {T <: AbstractFloat}
    Emn = energyspectrum(lx,ly,nx,ny,u,Λ=Λ)
    Vxy = vorticity(nx,ny,u,Λ=Λ)
    Uxy = zonalvelocity(lx,ly,nx,ny,u,Λ=Λ)
    Vyt = zonalvorticity(nx,ny,u)
    Dict("Emn"=>Emn,"Vxy"=>Vxy,"Uxy"=>Uxy,"Vyt"=>Vyt)
end

dumpstats(prob,u::Array{DSSField{T},1}) where {T<:AbstractFloat} = Dict("mEVs"=>modalevs(prob,u))

function dumpadjacency(lx::T,ly::T,nx::Int,ny::Int;fs::String,Λ::Int=nx-1) where {T <: AbstractFloat}
    A,C = adjacency(lx,ly,nx,ny,Λ=Λ)
    d = Dict("A"=>A,"C"=>C)
    NPZ.npzwrite(fs*".npz",d)
end
