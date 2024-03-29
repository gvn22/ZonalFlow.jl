function Base.write(prob,eqs,sol;dn::String="./",fn::String)
    mkpath(dn)
    @info "saving to" dn*fn*".jld2" 
    @save dn*fn*".jld2" sol
    NPZ.npzwrite(dn*fn*".npz",merge(dumpscalars(prob,sol),
                                    dumpfields(prob,sol),
                                    dumpstats(prob,eqs,sol),
                                    dumpcoeffs(prob,eqs,sol)
                                    ))
end

Base.write(prob,eqs::Vector{AbstractEquations},sols;dn::String="./",labels::Vector{String}=label(eqs)) = foreach(x->write(prob,x[1],x[2],dn=dn,fn=x[3]),zip(eqs,sols,labels))

tonpz(u) = reshape(cat(u...,dims=length(size(u[1]))),size(u[1])...,length(u))

function dumpcoeffs(prob,eqs,sol)
    x = zeros(eqs,prob.d)
    x .= fcoeffs(prob,eqs)
    F = forcingspectrum(prob.d,x)
    Dict("t"=>sol.t,"F"=>F)
end

function dumpnoise(prob,eqs,sol)
    (nx,ny) = size(prob.d)
    x = zeros(eqs,prob.d)
    x .= fcoeffs(prob,eqs)
    for m1=1:nx-1
        for n1 = -ny+1:ny-1
            x[n1+ny,m1+1] = x[n1+ny,m1+1]*sol.W[end][n1+ny,m1+1]
        end
    end
    W = resolvedfield(prob.d,x)
    Dict("t"=>sol.t,"W"=>W)
end

function dumpscalars(prob,sol;t0=18000.0)
    t0 = min(sol.t[end]/2.0,t0)
    Et = energy.(Ref(prob.d),sol.u)
    Zt = enstrophy.(Ref(prob.d),sol.u)
    Emt = zonalenergy.(Ref(prob.d),sol.u) |> tonpz
    Emtav = zonalenergy.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    Dict("t"=>sol.t,"Et"=>Et,"Emt"=>Emt,"Emtav"=>Emtav)
end

function dumpfields(prob,sol;t0=18000.0)
    t0 = min(sol.t[end]/2.0,t0)
    Emn = energyspectrum.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    Vxyav = vorticity.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    Vxy = vorticity.(Ref(prob.d),sol.u) |> tonpz
    Uxy = xvelocity.(Ref(prob.d),sol.u) |> tonpz
    Vyt = zonalvorticity.(Ref(prob.d),sol.u) |> tonpz
    Vytav = zonalvorticity.(Ref(prob.d),sol.u) |> x->timeaverage(sol.t,x,t0=t0) |> tonpz
    Dict("t"=>sol.t,"Emn"=>Emn,"Vxy"=>Vxy,"Vxyav"=>Vxyav,"Uxy"=>Uxy,"Vyt"=>Vyt,"Vytav"=>Vytav)
end

dumpstats(prob,eqs,sol) = Dict("mEVs"=> modaleigvals.(Ref(prob.d),sol.u) |> tonpz)

function dumpadjacency(prob,eqs::Union{NL,GQL};fn::String)
    A,B,C = adjacency(prob,eqs)
    NPZ.npzwrite(fn*".npz",Dict("A"=>A,"B"=>B,"C"=>C))
end

function dumpadjacency(prob,eqs::Union{NL,GQL},u;fn::String)
    A,B,C,Cl,Ch = adjacency(prob,eqs,u)
    NPZ.npzwrite(fn*".npz",Dict("A"=>A,"B"=>B,"C"=>C,"Cl"=>Cl,"Ch"=>Ch))
end
