function dumpenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u;fs::String,t0::Float64=200.0)

    Et,Zt = energy(lx,ly,nx,ny,u)
    Etav,Ztav = energy(lx,ly,nx,ny,t,u,t0=t0)

    Emt,Zmt = zonalenergy(lx,ly,nx,ny,u)
    Emtav,Zmtav = zonalenergy(lx,ly,nx,ny,t,u,t0=t0)

    d = Dict("t"=>t,"Et"=>Et,"Etav"=>Etav,"Emt"=>Emt,"Emtav"=>Emtav)
    NPZ.npzwrite(fs*".npz",d)
    nothing

end

function dumpfields(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u::Array{DNSField{T},1};fs::String,Λ::Int=nx-1) where T

    E = fourierenergy(lx,ly,nx,ny,u)
    E2 = energyspectrum(lx,ly,nx,ny,u,Λ=Λ)

    V = inversefourier(nx,ny,u)
    V2 = vorticity(nx,ny,u,Λ=Λ)

    U = zonalvelocity(lx,ly,nx,ny,u)
    U2 = zonalvelocity(lx,ly,nx,ny,u,Λ=Λ)
    #
    # W = zonalvorticity(lx,ly,nx,ny,u)
    # W2 = meanvorticity(lx,ly,nx,ny,u,Λ=Λ)

    # u1 = velocity(lx,ly,nx,ny,u)
    # u2 = inversefourier(nx,ny,u1)

    d = Dict("t"=>t,"E"=>E,"E2"=>E2,"V"=>V,"V2"=>V2)
    NPZ.npzwrite(fs*".npz",d)

    nothing

end

function dumpfields(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u::Array{GSSField{T},1};fs::String,Λ::Int) where T

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

function dumpstats()
end
