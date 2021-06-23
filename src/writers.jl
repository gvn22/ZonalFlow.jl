function Base.write(prob,eqs,sol;dn::String="",fn::String)
    mkpath(dn)
    NPZ.npzwrite(dn*fn*".npz",merge(dumpscalars(prob,sol),
                                    dumpfields(prob,sol),
                                    dumpstats(prob,eqs,sol)))
end

Base.write(prob,eqs::Vector{AbstractEquations},sols;dn::String,labels::Vector{String}=label(eqs)) = foreach(x->write(prob,x[1],x[2],dn=dn,fn=x[3]),zip(eqs,sols,labels))

tonpz(u) = reshape(cat(u...,dims=length(size(u[1]))),size(u[1])...,length(u))

function dumpscalars(prob,sol;t0=200.0)
    Et = energy.(Ref(prob.d),sol.u)
    Zt = enstrophy.(Ref(prob.d),sol.u)
    Emt = zonalenergy.(Ref(prob.d),sol.u) |> tonpz
    Emtav = zonalenergy.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x) |> tonpz
    Dict("t"=>sol.t,"Et"=>Et,"Emt"=>Emt,"Zt"=>Zt,"Emtav"=>Emtav)
end

function dumpfields(prob,sol)
    Emn = energyspectrum.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x) |> tonpz
    Vxy = vorticity.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x) |> tonpz
    Uxy = xvelocity.(Ref(prob.d),sol.u) |> tonpz
    Vyt = zonalvorticity.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x) |> tonpz
    Dict("t"=>sol.t,"Emn"=>Emn,"Vxy"=>Vxy,"Uxy"=>Uxy,"Vyt"=>Vyt)
end

dumpstats(prob,eqs,sol) = Dict("empty"=>0)
dumpstats(prob,eqs::CE2,sol) = Dict("mEVs"=> modaleigvals.(Ref(prob.d),sol.u) |> tonpz)

function dumpadjacency(lx::T,ly::T,nx::Int,ny::Int;fs::String,Λ::Int=nx-1) where {T <: AbstractFloat}
    A,C = adjacency(lx,ly,nx,ny,Λ=Λ)
    d = Dict("A"=>A,"C"=>C)
    NPZ.npzwrite(fs*".npz",d)
end
