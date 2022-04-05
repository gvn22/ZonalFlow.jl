"""
    Solve method for DiffEq solution
"""

function get_de_ic(prob,eqs,u0=nothing)
    if u0 == nothing
        @info "Setting a random initial condition..."
        # @info "Setting higher power IC"
        # rand(eqs,prob.d)
        # rand(eqs,prob.d,1e-5)
        # u0 = rand(eqs,prob.d,1e-4)
        # small noise to m = 1
        u0 = rand(eqs,prob.d,1e-4)
        # @info "with m = 3 set to 0..."
        # # u0 = convert(eqs,u0,prob.d)
        # for n = -prob.d.ny+1:prob.d.ny-1
        #     m = 3
        #     # u0[n+prob.d.ny,m+1] = 0.0
        #     u0.x[2][n+prob.d.ny,n+prob.d.ny,m] = 0.0
        # end
        return u0
    elseif typeof(u0) <: Number
        typeof(eqs) == CE2 ? fullrank(eqs,prob.d,u0) : rand(eqs,prob.d,u0)
    else
        @info "Converting QL solution to CE2 initial condition..."
        u0 = convert(eqs,u0,prob.d)
        # @info "Setting m = 3 to 0"
        # for n = -prob.d.ny+1:prob.d.ny-1
        #     m = 3
        #     u0.x[2][n+prob.d.ny,n+prob.d.ny,m] = 0.0
        # end
        # # small noise to m = 1
        # u0 = convert(eqs,u0,prob.d)
        # @info "Perturbing m = 2 in QL initial condition..."
        # for n = -prob.d.ny+1:prob.d.ny-1
        #     m = 2
        #     u0[n+prob.d.ny,m+1] = 1e-4*exp(im*rand(Uniform(0,2π)))
        #     for m = 3:prob.d.nx-1
        #         u0[n+prob.d.ny,m+1] = 0.0
        #     end
        # end
        # for n = 1:2prob.d.ny-1
        #     for n2 = 1:2prob.d.ny-1
        #         for m = 1:prob.d.nx-1
        #             u0.x[2][n2,n,m] += 1e-12*exp(im*rand(Uniform(0,2π)))
        #         end
        #         u0.x[2][n2,n,2] = 0.0
        #     end
        # end
        return u0
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
        poscondition(u,t,integrator) = t ∈ poschecktimes
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
