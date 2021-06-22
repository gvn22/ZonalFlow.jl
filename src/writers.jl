function Base.write(prob,eqs,sol;dn::String="",fn::String)
    mkpath(dn)
    NPZ.npzwrite(dn*fn*".npz",merge(dumpscalars(prob,eqs,sol),
                                    dumpfields(prob,eqs,sol),
                                    dumpstats(prob,eqs,sol)))
end

Base.write(prob,eqs::Vector{AbstractEquations},sols;dn::String,labels::Vector{String}=label(eqs)) = foreach(x->write(prob,x[1],x[2],dn=dn,fn=x[3]),zip(eqs,sols,labels))

tonpz(u) = reshape(cat(u...,dims=length(size(u[1]))),size(u[1])...,length(u))

function dumpscalars(prob,eqs,sol;t0=200.0)
    (lx,ly),(nx,ny) = length(prob.d),size(prob.d)
    Et,Zt = energy(lx,ly,nx,ny,sol.u)
    Etav,Ztav = energy(lx,ly,nx,ny,sol.t,sol.u,t0=t0)
    Emt,Zmt = zonalenergy(lx,ly,nx,ny,sol.u)
    Emtav,Zmtav = zonalenergy(lx,ly,nx,ny,sol.t,sol.u,t0=t0)
    Dict("t"=>sol.t,"Zt"=>Zt,"Ztav"=>Ztav,"Et"=>Et,"Etav"=>Etav,"Emt"=>Emt,"Emtav"=>Emtav)
end

function dumpfields(prob,eqs,sol)
    Emn = energyspectrum.(Ref(prob.d),sol.u) |> tonpz
    Vxy = vorticity.(Ref(prob.d),sol.u) |> tonpz
    Uxy = xvelocity.(Ref(prob.d),sol.u) |> tonpz
    Vyt = zonalvorticity.(Ref(prob.d),sol.u) |> tonpz
    Dict("t"=>sol.t,"Emn"=>Emn,"Vxy"=>Vxy,"Uxy"=>Uxy,"Vyt"=>Vyt)
end

dumpstats(prob,eqs,sol) = Dict("empty"=>0)
dumpstats(prob,eqs::CE2,sol) = Dict("mEVs"=>modaleigvals(prob.d,sol.u))

function dumpadjacency(lx::T,ly::T,nx::Int,ny::Int;fs::String,Λ::Int=nx-1) where {T <: AbstractFloat}
    A,C = adjacency(lx,ly,nx,ny,Λ=Λ)
    d = Dict("A"=>A,"C"=>C)
    NPZ.npzwrite(fs*".npz",d)
end
