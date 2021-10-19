function Base.write(prob,eqs,sol;dn::String="",fn::String)
    mkpath(dn)
    NPZ.npzwrite(dn*fn*".npz",merge(dumpscalars(prob,sol),
                                    dumpfields(prob,sol),
                                    dumpstats(prob,eqs,sol),
                                    dumpcoeffs(prob,eqs,sol),
                                    ))
end

Base.write(prob,eqs::Vector{AbstractEquations},sols;dn::String,labels::Vector{String}=label(eqs)) = foreach(x->write(prob,x[1],x[2],dn=dn,fn=x[3]),zip(eqs,sols,labels))

tonpz(u) = reshape(cat(u...,dims=length(size(u[1]))),size(u[1])...,length(u))

function dumpcoeffs(prob,eqs::CE2,sol)
    x = zeros(eqs,prob.d)
    x.x[2] .= fcoeffs(prob,eqs)
    F = forcingspectrum(prob.d,x)
    Dict("t"=>sol.t,"F"=>F)
end

function dumpcoeffs(prob,eqs,sol)
    x = zeros(eqs,prob.d)
    x .= fcoeffs(prob,eqs)
    F = forcingspectrum(prob.d,x)
    Dict("t"=>sol.t,"F"=>F)
end


function dumpscalars(prob,sol;t0=500.0)
    Et = energy.(Ref(prob.d),sol.u)
    Zt = enstrophy.(Ref(prob.d),sol.u)
    Emt = zonalenergy.(Ref(prob.d),sol.u) |> tonpz
    Emtav = zonalenergy.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    Dict("t"=>sol.t,"Et"=>Et,"Emt"=>Emt,"Emtav"=>Emtav,"Zt"=>Zt)
    # ! study eigenvalue structures for each zonal mode in second cumulant
    # Em0nt = modalenergy.(Ref(prob.d),sol.u,m=0) |> tonpz
    # Em1nt = modalenergy.(Ref(prob.d),sol.u,m=1) |> tonpz
    # Em2nt = modalenergy.(Ref(prob.d),sol.u,m=2) |> tonpz
    # Em0ntav = modalenergy.(Ref(prob.d),sol.u,m=0) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    # Em1ntav = modalenergy.(Ref(prob.d),sol.u,m=1) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    # Em2ntav = modalenergy.(Ref(prob.d),sol.u,m=2) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    # Dict("t"=>sol.t,"Et"=>Et,"Emt"=>Emt,"Zt"=>Zt,"Emtav"=>Emtav,"Em0nt"=>Em0nt,"Em1nt"=>Em1nt,"Em2nt"=>Em2nt,
    #     "Em0ntav"=>Em0ntav,"Em1ntav"=>Em1ntav,"Em2ntav"=>Em2ntav)
end

function dumpfields(prob,sol;t0=500.0)
    Emn = energyspectrum.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    Fyym1 = secondcumulant.(Ref(prob.d),sol.u,m=1) |> tonpz
    Fyym2 = secondcumulant.(Ref(prob.d),sol.u,m=2) |> tonpz
    Fyym3 = secondcumulant.(Ref(prob.d),sol.u,m=3) |> tonpz
    Vxyav = vorticity.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    Vxy = vorticity.(Ref(prob.d),sol.u) |> tonpz
    Uxy = xvelocity.(Ref(prob.d),sol.u) |> tonpz
    Vyt = zonalvorticity.(Ref(prob.d),sol.u) |> tonpz
    Vytav = zonalvorticity.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    Dict("t"=>sol.t,"Emn"=>Emn,"Vxy"=>Vxy,"Vxyav"=>Vxyav,"Uxy"=>Uxy,"Vyt"=>Vyt,"Vytav"=>Vytav,"Fyym1"=>Fyym1,"Fyym2"=>Fyym2,"Fyym3"=>Fyym3)
end

dumpstats(prob,eqs,sol) = Dict("empty"=>0)
dumpstats(prob,eqs::GQL,sol) = Dict("mEVs"=> convert.(Ref(CE2()),sol.u,Ref(prob.d)) |> x-> modaleigvals.(Ref(prob.d),x) |> tonpz)
dumpstats(prob,eqs::CE2,sol) = Dict("mEVs"=> modaleigvals.(Ref(prob.d),sol.u) |> tonpz)
# !!! time averaged second cumulants from QL/CE2
# dumpstats(prob,eqs::GQL,sol) = Dict("mEVs"=> convert.(Ref(CE2()),sol.u,Ref(prob.d)) |> x-> timeaverage(sol.t,x,t0=500.0) |> x-> modaleigvals.(Ref(prob.d),x) |> tonpz)
# dumpstats(prob,eqs::CE2,sol) = Dict("mEVs"=> timeaverage(sol.t,sol.u,t0=500.0) |> x-> modaleigvals.(Ref(prob.d),x) |> tonpz)

function dumpadjacency(prob,eqs::Union{NL,GQL};fn::String)
    A,B,C = adjacency(prob,eqs)
    NPZ.npzwrite(fn*".npz",Dict("A"=>A,"B"=>B,"C"=>C))
end

function dumpadjacency(prob,eqs::Union{NL,GQL},sol;fn::String)
    A,B,C = adjacency(prob,eqs,sol)
    NPZ.npzwrite(fn*".npz",Dict("A"=>A,"B"=>B,"C"=>C))
end
