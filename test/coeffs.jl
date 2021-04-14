using Revise
using Test
using ZonalFlow
using BenchmarkTools

lx = 2.0π
ly = 2.0π
nx = rand(collect(2:64))
ny = rand(collect(2:64))

Ω = 2.0π
θ = rand(0.0:90.0)
β = 2*Ω*cos(deg2rad(θ))
μ = rand(1.0:100.0)*0.01
ν = rand(1.0:100.0)*0.01
ν₄ = rand(1.0:100.0)*0.01

dom = Domain(extent=(lx,ly),res=(nx,ny))
coeffs = Coefficients(Ω=Ω,θ=θ,μ=μ,ν=ν,ν₄=ν₄)

Ξ = rand(1.0:10.0)*0.1
Δθ = rand(1.0:5.0)*0.04
τ = 1.0/μ

forcing = PointJet(Ξ=Ξ,Δθ=Δθ,τ=τ)
prob = BetaPlane(dom,coeffs,forcing)

@benchmark acoeffs(prob)
@benchmark acoeffs(ly,ny,Ξ*Ω,Δθ,τ)
A = acoeffs(prob)
A0 = acoeffs(ly,ny,Ξ*Ω,Δθ,τ)
@test isapprox(A,A0)
@code_warntype acoeffs(prob)

@benchmark bcoeffs(prob.d,prob.c)
@benchmark bcoeffs(length(prob.d)...,size(prob.d)...,β,μ,ν,ν₄)
# @code_warntype bcoeffs(prob.d,prob.c)
B = bcoeffs(prob.d,prob.c)
B0 = bcoeffs(length(prob.d)...,size(prob.d)...,β,μ,ν,ν₄)
@test isapprox(B,B0)

@benchmark ccoeffs(prob.d,NL())
@benchmark ccoeffs(lx,ly,nx,ny)
# @code_warntype ccoeffs(prob.d,NL())
Cp,Cm = ccoeffs(prob.d,NL())
Cp0,Cm0 = ccoeffs(lx,ly,nx,ny)
@test isapprox(Cp,Cp0)
@test isapprox(Cm,Cm0)

Λ = rand(0:nx-1)
Cp,Cm = ccoeffs(prob.d,GQL(Λ=Λ))
Cp0,Cm0 = ccoeffs(lx,ly,nx,ny,Λ)
@test isapprox(Cp,Cp0)
@test isapprox(Cm,Cm0)
