"""
    Solve method for DiffEq solution
"""
function get_de_ic(prob,eqs,u0=nothing)
    if u0 == nothing
        @info "Setting a random initial condition..."
        rand(eqs,prob.d)
    elseif typeof(u0) <: Number
        @info "Setting a DNS initial condition..."
        rand(eqs,prob.d,u0)
    else
        @info "Converting QL solution to CE2 initial condition..."
        convert(eqs,u0,prob.d)
    end
end

# function get_de_ic(prob::BetaPlane{T,Kolmogorov{T}},eqs::Union{NL,GQL},u0=nothing) where T
#     a = rand(eqs,prob.d)
#     a[:,1] += acoeffs(prob)
#     u0 == nothing ? a : convert(eqs,u0,prob.d)
# end

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

get_de_probalg(prob,eqs,u0,t,p) = ODEProblem(f!,u0,t,p), RK4()
get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs::CE2,u0,t,p) where T = ODEProblem(f!,u0,t,p), Heun()
get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs,u0,t,p) where T = SDEProblem(f!,g!,u0,t,p), SRIW1()

get_de_kwargs(prob,eqs::AbstractEquations,tspan;kwargs...) = kwargs

# function get_de_kwargs(prob::BetaPlane{T,Kolmogorov{T}},eqs::CE2,tspan;kwargs...) where T<:AbstractFloat
#     info = SavedValues(T,Tuple{Int,Array{T,1},Array{Int,1},Array{T,2}})
#     func(u,t,integrator) = rankis(integrator.u.x[2],prob.d.nx,prob.d.ny)
#     rankchecktimes = [i for i=0.0:10.0:tspan[2]]
#     rankcheckcallback = SavingCallback(func,info,save_start=true,save_everystep=false,saveat=rankchecktimes)
#     merge((callback=rankcheckcallback,tstops=rankchecktimes),kwargs),info
# end

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
    if(!eqs.poscheck) return kwargs end
    @info "Setting positivity check callback..."
    poschecktimes = [tt for tt=tspan[1]:eqs.poscheckat:tspan[2]]
    temp = zeros(ComplexF64,2*prob.d.ny-1,prob.d.nx-1,2*prob.d.ny-1,prob.d.nx-1)
    condition(u,t,integrator) = t ∈ poschecktimes
    affect!(integrator) = positivity!(integrator.u.x[2],temp,prob.d.nx,prob.d.ny)
    poscheckcb = DiscreteCallback(condition,affect!,save_positions=(false,false))
    merge((callback=poscheckcb,tstops=poschecktimes),kwargs)
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
