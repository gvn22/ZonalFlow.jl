function dumpenergy(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u;fs::String,t0::Float64=200.0)

    Et,Zt = energy(lx,ly,nx,ny,u)
    Etav,Ztav = energy(lx,ly,nx,ny,t,u,t0=t0)

    Emt,Zmt = zonalenergy(lx,ly,nx,ny,u)
    Emtav,Zmtav = zonalenergy(lx,ly,nx,ny,t,u,t0=t0)

    d = Dict("t"=>t,"Et"=>Et,"Etav"=>Etav,"Emt"=>Emt,"Emtav"=>Emtav)
    NPZ.npzwrite(fs*".npz",d)
    nothing

end

function dumpfields(lx::Float64,ly::Float64,nx::Int,ny::Int,t::Array{Float64,1},u;fs::String)

    E = fourierenergy(lx,ly,nx,ny,u)
    V = inversefourier(nx,ny,u)
    U = zonalvelocity(lx,ly,nx,ny,u)

    d = Dict("t"=>t,"E"=>E,"V"=>V,"U"=>U)
    NPZ.npzwrite(fs*".npz",d)

    nothing

end

function dumpstats()
end
