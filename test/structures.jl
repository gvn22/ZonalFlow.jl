using ZonalFlow
using Test

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
Δ = 0.1
forcing = Stochastic(kf=kf,dk=dk,ε=ε,Δ = 0.1)
@test forcing.kf == kf
@test forcing.dk == dk
@test forcing.ε == ε

prob1 = BetaPlane(dom,coeffs,forcing)
prob2 = BetaPlane(coeffs,forcing,extent=(lx,ly),res=(nx,ny))
@test prob1.d == dom
@test prob1.f == forcing
@test prob1.c == coeffs
@test prob2 == prob1

@test zeros(GQL(rand(0:dom.nx-1)),dom) == zeros(NL(),dom)

Λ = 2

A = zeros(GCE2(Λ),dom)
@test isequal(A.l,A.x[1])
@test isequal(A.h,A.x[2])

B = zeros(CE2(),dom)
@test isequal(B.l,B.x[1])
@test isequal(B.h,B.x[2])
