using ZonalFlow
using Test

domain  = Domain(extent=(2π,2π),res=(4,4));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.0,ν=0.0,ν₄=1.0,linear=true);
forcing = Stochastic(kf=3,dk=1,ε=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=500
           );

eqs     = [NL(),CE2()];
eqs     = append!(eqs,[GQL(Λ=l) for l=0:prob.d.nx-1])
eqs     = append!(eqs,[GCE2(Λ=l) for l=0:prob.d.nx-1])
sols    = integrate(prob,eqs,tspan;tsargs...);

@testset "Linear Equations" begin
        for sol in sols
                E,Z = energy(length(prob.d)...,size(prob.d)...,sol.t,sol.u);
                @test E[1] == E[end]
                @test Z[1] == Z[end]
        end
end

coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.0,ν=0.0,ν₄=1.0,linear=false);
prob    = BetaPlane(domain,coeffs,forcing);
sols    = integrate(prob,eqs,tspan;tsargs...);

@testset "Nonlinear Equations" begin
        for sol in sols
                E,Z = energy(length(prob.d)...,size(prob.d)...,sol.t,sol.u);
                @test E[1] == E[end]
                @test Z[1] == Z[end]
        end
end
