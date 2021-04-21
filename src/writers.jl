function Base.write(prob,eqs::Vector{AbstractEquations},sols;dn::String,labels::Vector{String}=label(eqs))
    foreach(x->write(prob,x[1],x[2],dn=dn,fn=x[3]),zip(eqs,sols,labels))
end

function Base.write(prob,eqs,sol;dn::String,fn::String)
    (lx,ly),(nx,ny),Λ = length(prob.d),size(prob.d),lambda(prob,eqs)
    t,u = sol.t,sol.u

    Et,Zt = energy(lx,ly,nx,ny,u)
    Etav,Ztav = energy(lx,ly,nx,ny,t,u)
    Emt,Zmt = zonalenergy(lx,ly,nx,ny,u)
    Emtav,Zmtav = zonalenergy(lx,ly,nx,ny,t,u)

    # Emn = fourierenergy(lx,ly,nx,ny,u)
    # Vxy = inversefourier(nx,ny,u)
    # Emn2 = energyspectrum(lx,ly,nx,ny,u,Λ=Λ)
    # Vxy2 = vorticity(nx,ny,u,Λ=Λ)
    # Uxy = zonalvelocity(lx,ly,nx,ny,u)
    # Uxy2 = zonalvelocity(lx,ly,nx,ny,u,Λ=Λ)

    d = Dict("t"=>t,"Zt"=>Zt,"Ztav"=>Ztav,"Et"=>Et,"Etav"=>Etav,"Emt"=>Emt,"Emtav"=>Emtav)
    mkpath(dn)
    NPZ.npzwrite(dn*fn*".npz",d)
    nothing
end

function dumpenergy(lx::T,ly::T,nx::Int,ny::Int,t::Array{T,1},u;fs::String,t0::Float64=200.0) where {T <: AbstractFloat}

    Et,Zt = energy(lx,ly,nx,ny,u)
    Etav,Ztav = energy(lx,ly,nx,ny,t,u,t0=t0)

    Emt,Zmt = zonalenergy(lx,ly,nx,ny,u)
    Emtav,Zmtav = zonalenergy(lx,ly,nx,ny,t,u,t0=t0)

    d = Dict("t"=>t,"Et"=>Et,"Etav"=>Etav,"Emt"=>Emt,"Emtav"=>Emtav)
    NPZ.npzwrite(fs*".npz",d)
    nothing

end

function dumpfields(lx::T,ly::T,nx::Int,ny::Int,t::Array{Float64,1},u::Array{DNSField{T},1};fs::String,Λ::Int=nx-1) where {T <: AbstractFloat}

    E = fourierenergy(lx,ly,nx,ny,u)
    E2 = energyspectrum(lx,ly,nx,ny,u,Λ=Λ)

    V = inversefourier(nx,ny,u)
    V2 = vorticity(nx,ny,u,Λ=Λ)

    U = zonalvelocity(lx,ly,nx,ny,u)
    U2 = zonalvelocity(lx,ly,nx,ny,u,Λ=Λ)

    # W = zonalvorticity(lx,ly,nx,ny,u)
    # W2 = meanvorticity(lx,ly,nx,ny,u,Λ=Λ)

    # u1 = velocity(lx,ly,nx,ny,u)
    # u2 = inversefourier(nx,ny,u1)

    d = Dict("t"=>t,"E"=>E,"E2"=>E2,"V"=>V,"V2"=>V2)
    NPZ.npzwrite(fs*".npz",d)

    nothing

end

function dumpfields(lx::T,ly::T,nx::Int,ny::Int,t::Array{Float64,1},u::Array{GSSField{T},1};fs::String,Λ::Int) where {T <: AbstractFloat}

    E = fourierenergy(lx,ly,nx,ny,Λ,u)
    E2 = energyspectrum(lx,ly,nx,ny,u,Λ=Λ)

    V = inversefourier(nx,ny,Λ,u)
    V2 = vorticity(nx,ny,u,Λ=Λ)
    # U = zonalvelocity(lx,ly,nx,ny,u)

    # u1 = velocity(lx,ly,nx,ny,u)
    # u2 = inversefourier(nx,ny,u1)

    d = Dict("t"=>t,"E"=>E,"E2"=>E2,"V"=>V,"V2"=>V2)
    NPZ.npzwrite(fs*".npz",d)

    nothing

end

function dumpstats() end

function dumpadjacency(lx::T,ly::T,nx::Int,ny::Int;fs::String,Λ::Int=nx-1) where {T <: AbstractFloat}
    A,C = adjacency(lx,ly,nx,ny,Λ=Λ)
    d = Dict("A"=>A,"C"=>C)
    NPZ.npzwrite(fs*".npz",d)
    nothing
end
