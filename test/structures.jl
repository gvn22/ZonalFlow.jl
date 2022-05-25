using Test
using ZonalFlow
using BenchmarkTools

lx = 2.0π
ly = 2.0π
nx = 12
ny = 12

dom = Domain(extent=(lx,ly),res=(nx,ny))
@test dom.lx == lx
@test dom.ly == ly
@test dom.nx == nx
@test dom.ny == ny

Ω = 2.0π
θ = 0.0
μ = 0.01
ν = 0.0
ν₄ = 1.0

coeffs = Coefficients(Ω=Ω,θ=θ,μ=μ,ν=ν,ν₄=ν₄)
@test coeffs.Ω == Ω
@test coeffs.θ == θ
@test coeffs.μ == μ
@test coeffs.ν == ν
@test coeffs.ν₄ == ν₄

Ξ = 1.0
Δθ = 0.05
τ = 10.0
forcing = PointJet(Ξ=Ξ,Δθ=Δθ,τ=τ)
@test forcing.Ξ == Ξ
@test forcing.Δθ == Δθ
@test forcing.τ == τ

A₁ = 1.0
A₄ = 2.0
forcing = Kolmogorov(A₁=A₁,A₄=A₄)
@test forcing.A₁ == A₁
@test forcing.A₄ == A₄

kf = 6
dk = 1
ε = 1.0e-3
forcing = Stochastic(kf=kf,dk=dk,ε=ε)
@test forcing.kf == kf
@test forcing.dk == dk
@test forcing.ε == ε

prob1 = BetaPlane(dom,coeffs,forcing)
prob2 = BetaPlane(coeffs,forcing,extent=(lx,ly),res=(nx,ny))
@test prob1.d == dom
@test prob1.f == forcing
@test prob1.c == coeffs
@test prob2 == prob1

@benchmark zeros(NL(),dom)
zeros(NL(),dom)
@test zeros(GQL(rand(0:dom.nx-1)),dom) == zeros(NL(),dom)

@benchmark zeros(CE2(),dom)
A = zeros(GCE2(2),dom)
@test isequal(A.l,A.x[1])
@test isequal(A.h,A.x[2])

B = zeros(CE2(),dom)
@test isequal(B.l,B.x[1])
@test isequal(B.h,B.x[2])

C = zeros(GCE2(2),dom)
@test isequal(C.l,C.x[1])
@test isequal(C.h,C.x[2])
