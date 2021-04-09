"""
    Solve method for DiffEq solution
"""
Base.zeros(::Field{T},nx,ny) where T<:AbstractFloat = zeros(Complex{T},2ny-1,nx)
Base.zeros(::FirstCumulant{T},nx,ny) where T<:AbstractFloat = zeros(Complex{T},2ny-1)
Base.zeros(::SecondCumulant{T},nx,ny) where T<:AbstractFloat = zeros(Complex{T},2ny-1,2ny-1,nx-1)
Base.zeros(::FieldBilinear{T},nx,ny,Λ) where T<:AbstractFloat = zeros(Complex{T},2ny-1,nx-Λ,2ny-1,nx-Λ)
# Random.rand(::Field{T},dims) = rand(Field,dims)...

get_de_ic(d::Domain{T},eqs) where T<:AbstractFloat = zeros(Field{T},size(d)...)
get_de_ic(d::Domain{T},eqs::CE2) where T<:AbstractFloat = ArrayPartition(zeros(::FirstCumulant{T},size(d)...),zeros(::SecondCumulant{T},size(d)...))
get_de_ic(d::Domain{T},eqs::GCE2) where T<:AbstractFloat = ArrayPartition(zeros(::Field{T},eqs.Λ+1,d.ny),zeros(::FieldBilinear{T},size(d)...,eqs.Λ)

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
    ODEProblem(f,u0,t,p), RK4()
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
    u0 = get_de_ic(prob.d,eqs)
    p = get_de_params(prob,eqs)
    _prob,_alg = get_de_probalg(prob,eqs,u0,p,tspan)
    solve(_prob,_alg,kwargs...)
end
