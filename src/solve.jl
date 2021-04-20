"""
    Solve method for DiffEq solution
"""
get_de_ic(prob,eqs) = zeros(eqs,prob.d)

function get_de_params(prob,eqs)::AbstractParams
    A = acoeffs(prob)
    B = bcoeffs(prob)
    C⁺,C⁻ = ccoeffs(prob,eqs)
    F = fcoeffs(prob,eqs)
    get_de_p(prob.d,eqs,[A,B,C⁺,C⁻,F])
end

get_de_p(d,eqs::NL,p) = NLParams(d.nx,d.ny,p...)
get_de_p(d,eqs::GQL,p) = GQLParams(d.nx,d.ny,eqs.Λ,p...)
get_de_p(d::Domain{T},eqs::CE2,p) where T = [d.nx,d.ny,p...,similar(FirstCumulant{T},d),similar(SecondCumulant{T},d),similar(SecondCumulant{T},d)]
get_de_p(d::Domain{T},eqs::GCE2,p) where T = [d.nx,d.ny,eqs.Λ,p...,similar(Field{T},d,Λ=eqs.Λ),similar(FieldBilinear{T},d,Λ=eqs.Λ),similar(FieldBilinear{T},d,Λ=eqs.Λ)]

get_de_eqs(::NL) = f!,g!
get_de_eqs(::GQL) = f!,g!
get_de_eqs(::GCE2) = gce2_eqs!,g!
get_de_eqs(::CE2) = ce2_eqs!

function get_de_probalg(prob,eqs,u0,t,p)
    f,g = get_de_eqs(eqs)
    ODEProblem(f,u0,t,p), RK4()
end

function get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs,u0,t,p) where T
    f,g = get_de_eqs(eqs)
    W0 = zeros(eqs,prob.d)
    Random.seed!(123)
    SDEProblem(f,g,u0,t,p,noise=noise!(t[1],W0)), EulerHeun()
end

function get_de_probalg(prob::BetaPlane{T,Stochastic{T}},eqs::CE2,u0,t,p) where T
    f = get_de_eqs(eqs)
    ODEProblem(f,u0,t,p), Heun()
end

function integrate(prob,eqs::AbstractEquations,tspan;kwargs...)
    u0 = get_de_ic(prob,eqs)
    p  = get_de_params(prob,eqs)
    _prob,_alg = get_de_probalg(prob,eqs,u0,tspan,p)
    @time solve(_prob,_alg;kwargs...)
end

integrate(prob,eqs::Vector{AbstractEquations},tspan;kwargs...) = [integrate(prob,eq,tspan;kwargs...) for eq in eqs]
