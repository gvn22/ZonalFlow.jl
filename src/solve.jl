"""
    Solve method for DiffEq solution
"""
function get_de_ic(prob,eqs,u0=nothing)
    if u0 == nothing
        @info "Setting a random initial condition..."
        rand(eqs,prob.d)
    elseif typeof(u0) <: Number
        typeof(eqs) == CE2 ? fullrank(eqs,prob.d,u0) : rand(eqs,prob.d,u0)
    else
        @info "Converting QL solution to CE2 initial condition..."
        convert(eqs,u0,prob.d)
        # u0 = convert(eqs,u0,prob.d)
        # small noise to m = 1
        # for n = 1:2prob.d.ny-1
        #     u0.x[2][n,n,1] += 1e-6*exp(im*rand(Uniform(0,2π)))
        #     @show u0.x[2][n,n,1]
        # end
        # return u0
    end
end

function get_de_params(prob,eqs)::AbstractParams
    A = acoeffs(prob)
    B = bcoeffs(prob)
    C⁺,C⁻ = ccoeffs(prob,eqs)
    F = fcoeffs(prob,eqs)
    get_de_p(prob.d,eqs,[A,B,C⁺,C⁻,F])
end

get_de_p(d,eqs::NL,p) = NLParams(d.nx,d.ny,p...)
get_de_p(d,eqs::GQL,p) = GQLParams(d.nx,d.ny,eqs.Λ,p...)
get_de_p(d,eqs::CE2,p) = CE2Params(d.nx,d.ny,p...)
get_de_p(d,eqs::GCE2,p) = GCE2Params(d.nx,d.ny,eqs.Λ,p...)

get_de_probalg(prob,eqs,u0,t,p) = ODEProblem(f!,u0,t,p), DP5()
get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs::CE2,u0,t,p) where T = ODEProblem(f!,u0,t,p), RK4()
get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs,u0,t,p) where T = SDEProblem(f!,g!,u0,t,p), SRIW1()

get_de_kwargs(prob,eqs::AbstractEquations,tspan;kwargs...) = kwargs

function get_de_kwargs(prob,eqs::GCE2,tspan;kwargs...)
    if(!eqs.poscheck) return kwargs end
    @info "Setting positivity check callback..."
    poschecktimes = [tt for tt=tspan[1]:eqs.poscheckat:tspan[2]]
    condition(u,t,integrator) = t ∈ poschecktimes && !ispositive(u.x[2],prob.d.nx,prob.d.ny,eqs.Λ)
    affect!(integrator) = positivity!(integrator.u.x[2],prob.d.nx,prob.d.ny,eqs.Λ)
    poscheckcb = DiscreteCallback(condition,affect!,save_positions=(false,false))
    merge((callback=poscheckcb,tstops=poschecktimes),kwargs)
end

function get_de_kwargs(prob,eqs::CE2,tspan;kwargs...)
    if(!eqs.poscheck && !eqs.eigmax)
        return kwargs
    elseif(eqs.poscheck && !eqs.eigmax)
        @info "Setting positivity check callback..."
        poschecktimes = [tt for tt=tspan[1]:eqs.poscheckat:tspan[2]]
        condition(u,t,integrator) = t ∈ poschecktimes
        temp = zeros(ComplexF64,2*prob.d.ny-1,prob.d.nx-1,2*prob.d.ny-1,prob.d.nx-1)
        posaffect!(integrator) = positivity!(integrator.u.x[2],temp,prob.d.nx,prob.d.ny)
        poscondition(u,t,integrator) = t > 0.0
        poscheckcb = DiscreteCallback(poscondition,posaffect!,save_positions=(false,false))
        return merge((callback=poscheckcb,tstops=poschecktimes),kwargs)
    else
        @info "Removing all but the largest modal eigenvalues..."
        temp = zeros(ComplexF64,2*prob.d.ny-1,2*prob.d.ny-1,prob.d.nx-1)
        maxcondition(u,t,integrator) = true
        maxaffect!(integrator) = truncatecumulant!(prob.d,integrator.u.x[2],temp)
        eigmaxcb = DiscreteCallback(maxcondition,maxaffect!,save_positions=(false,false))
        return merge((callback=eigmaxcb,tspan=tspan),kwargs)
    end
end

function integrate(prob,eqs::AbstractEquations,tspan;u0=nothing,kwargs...)
    Random.seed!(123)
    _u0 = get_de_ic(prob,eqs,u0)
    _p  = get_de_params(prob,eqs)
    _prob,_alg = get_de_probalg(prob,eqs,_u0,tspan,_p)
    _kwargs = get_de_kwargs(prob,eqs,tspan;kwargs...)
    @time solve(_prob,_alg;_kwargs...)
end

integrate(prob,eqs::Vector{AbstractEquations},tspan;kwargs...) = [integrate(prob,eq,tspan;kwargs...) for eq in eqs]
