"""
    Solve method for DiffEq solution
"""
# overload Base and Random methods for ICs
# Base.zeros(::NL,dims) = zeros(Field,dims)
# Random.rand(::NL,dims) = rand(Field,dims)...
# Base.show(::p) create a structure of type p

function get_de_params(prob,eqs)
    A = acoeffs(prob)
    B = bcoeffs(prob)
    C⁺,C⁻ = ccoeffs(prob,eqs)
    F = fcoeffs(prob,eqs)
    [A,B,C⁺,C⁻,F]
end

get_de_eqs(::NL) = nl_eqs!,unit_eqs!
get_de_eqs(::GQL) = gql_eqs!,unit_eqs!
get_de_eqs(::GCE2) = gce2_eqs!,unit_gce2_eqs!
get_de_eqs(::CE2) = ce2_eqs!

function get_de_probalg(prob,eqs,u0,p,t)
    f,g = get_de_eqs(eqs)
    ODEProblem(f,u0,p,t), RK4()
end

function get_de_probalg(prob::BetaPlane{Stochastic},eqs,u0,p,t)
    f... = get_de_eqs(eqs)
    W0 = zeros(eqs,size(prob.d))
    SDEProblem(f...,u0,t,p,noise=noise!(t[0],W0)), EulerHeun()
end

function get_de_probalg(prob::BetaPlane{Stochastic},eqs::CE2,u0,p,t)
    f = get_de_eqs(eqs)
    ODEProblem(f,u0,t,p), Heun()
end

function solve(prob,eqs;tspan,kwargs...)
    u0 = zeros(eqs,size(prob.d))
    p = get_de_params(prob,eqs)
    _prob,_alg = get_de_probalg(prob,eqs,u0,p,tspan)
    solve(_prob,_alg,kwargs...)
end
