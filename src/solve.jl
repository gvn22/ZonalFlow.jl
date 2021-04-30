"""
    Solve method for DiffEq solution
"""
function get_de_ic(prob,eqs,u0=nothing)
    Random.seed!(123)
    u0 == nothing ? rand(eqs,prob.d) : convert(eqs,u0,prob.d)
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

get_de_probalg(prob,eqs,u0,t,p) = ODEProblem(f!,u0,t,p), RK4()
get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs::CE2,u0,t,p) where T = ODEProblem(f!,u0,t,p), Heun()
get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs,u0,t,p) where T =  SDEProblem(f!,g!,u0,t,p), EulerHeun()
# function get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs,u0,t,p) where T
    # W0 = zeros(eqs,prob.d)
    # Random.seed!(123)
    # SDEProblem(f!,g!,u0,t,p,noise=noise!(t[1],W0)), EulerHeun()
     # SDEProblem(f!,g!,u0,t,p), EulerHeun()
# end

function integrate(prob,eqs::AbstractEquations,tspan;u0=nothing,kwargs...)
    u0 = get_de_ic(prob,eqs,u0)
    p  = get_de_params(prob,eqs)
    _prob,_alg = get_de_probalg(prob,eqs,u0,tspan,p)
    @time solve(_prob,_alg;kwargs...)
end

integrate(prob,eqs::Vector{AbstractEquations},tspan;kwargs...) = [integrate(prob,eq,tspan;kwargs...) for eq in eqs]
