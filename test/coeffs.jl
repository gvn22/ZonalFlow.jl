using Revise
using Test
using ZonalFlow

lx = 2.0π
ly = 2.0π
nx = rand(collect(2:64))
ny = rand(collect(2:64))
dom = Domain(extent=(lx,ly),res=(nx,ny))

Ω = 2.0π
θ = rand(0.0:90.0)
β = 2*Ω*cos(deg2rad(θ))
μ = rand(1.0:100.0)*0.01
ν = rand(1.0:100.0)*0.01
ν₄ = rand(1.0:100.0)*0.01
linear = [true,false][rand(1:2)]
coeffs = Coefficients(Ω=Ω,θ=θ,μ=μ,ν=ν,ν₄=ν₄,linear=linear)

Λ = rand(0:nx-1)

@testset "PointJet" begin

    Ξ = rand(1.0:10.0)*0.1
    Δθ = rand(1.0:5.0)*0.04
    τ = 1.0/μ
    forcing = PointJet(Ξ=Ξ,Δθ=Δθ,τ=τ)
    prob = BetaPlane(dom,coeffs,forcing)

    A = acoeffs(prob)
    A0 = acoeffs(ly,ny,Ξ*Ω,Δθ,τ)
    @test isapprox(A,A0)

    B = bcoeffs(prob)
    B0 = bcoeffs(length(prob.d)...,size(prob.d)...,β,μ,ν,ν₄)
    @test isapprox(B,B0)

    Cp,Cm = ccoeffs(prob,NL())
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,GQL(Λ))
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,Λ)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,GCE2(Λ))
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,Λ)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,CE2())
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,0)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    F = fcoeffs(prob,NL())
    F0 = fcoeffs(nx,ny)
    @test isequal(F,F0)

    F = fcoeffs(prob,GQL(Λ))
    F0 = fcoeffs(nx,ny)
    @test isequal(F,F0)

    F = fcoeffs(prob,GCE2(Λ))
    F0 = fcoeffs(nx,ny,Λ)
    @test isapprox(F,F0)

    F = fcoeffs(prob,CE2())
    F0 = fcoeffs(nx,ny,0)
    @test isequal(F,F0)

end



@testset "Kolmogorov" begin
    A₁ = -rand(1.0:10.0)
    A₄ = -rand(1.0:10.0)

    g = zeros(ComplexF64,2*ny-1);
    g[1+ny] = g[-1+ny]= A₁;
    g[4+ny] = g[-4+ny]= 4.0*A₄;

    coeffs = Coefficients(Ω=Ω,θ=θ,μ=μ,ν=ν,ν₄=ν₄,linear=linear)
    forcing = Kolmogorov(A₁=A₁,A₄=A₄)
    prob = BetaPlane(dom,coeffs,forcing)

    A = acoeffs(prob)
    A0 = acoeffs(ly,ny,g)
    @test isapprox(A,A0)

    B = bcoeffs(prob)
    B0 = bcoeffs(length(prob.d)...,size(prob.d)...,β,μ,ν,ν₄)
    @test isapprox(B,B0)

    Cp,Cm = ccoeffs(prob,NL())
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,GQL(Λ))
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,Λ)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,GCE2(Λ))
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,Λ)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,CE2())
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,0)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    F = fcoeffs(prob,NL())
    F0 = fcoeffs(nx,ny)
    @test isequal(F,F0)

    F = fcoeffs(prob,GQL(Λ))
    F0 = fcoeffs(nx,ny)
    @test isequal(F,F0)

    F = fcoeffs(prob,GCE2(Λ))
    F0 = fcoeffs(nx,ny,Λ)
    @test isapprox(F,F0)

    F = fcoeffs(prob,CE2())
    F0 = fcoeffs(nx,ny,0)
    @test isequal(F,F0)

end

@testset "Stochastic" begin

    kf = Int(rand(floor(nx/3):floor(2nx/3)))
    dk = rand(1:3)
    ε = rand(1.0:10.0)*1e-3

    forcing = Stochastic(kf=kf,dk=dk,ε=ε)
    prob = BetaPlane(dom,coeffs,forcing)

    A = acoeffs(prob)
    A0 = acoeffs(ny)
    @test isapprox(A,A0)

    B = bcoeffs(prob)
    B0 = bcoeffs(length(prob.d)...,size(prob.d)...,β,μ,ν,ν₄)
    @test isapprox(B,B0)

    Cp,Cm = ccoeffs(prob,NL())
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,GQL(Λ))
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,Λ)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,GCE2(Λ))
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,Λ)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    Cp,Cm = ccoeffs(prob,CE2())
    Cp0,Cm0 = linear == true ? ccoeffs(nx,ny) : ccoeffs(lx,ly,nx,ny,0)
    @test isapprox(Cp,Cp0)
    @test isapprox(Cm,Cm0)

    F = fcoeffs(prob,NL())
    F0 = fcoeffs(nx,ny,kf,dk,ε)
    @test isapprox(F,F0)

    F = fcoeffs(prob,GQL(Λ=Λ))
    F0 = fcoeffs(nx,ny,kf,dk,ε)
    @test isapprox(F,F0)

    F = fcoeffs(prob,GCE2(Λ=Λ))
    F0 = fcoeffs(nx,ny,Λ,kf,dk,ε)
    @test isapprox(F,F0)

    F = fcoeffs(prob,CE2())
    F0 = fcoeffs(nx,ny,0,kf,dk,ε)
    @test isapprox(F,F0)

end
