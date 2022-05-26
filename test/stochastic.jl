using ZonalFlow
using Test

tspan   = (0.0,200.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=false,
            progress_steps=10000,
            save_everystep=false,
            saveat=50
           );

domain  = Domain(extent=(2π,2π),res=(6,6));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0);
forcing = Stochastic(kf=3,dk=1,ε=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

eqs = [NL(),GQL(0),GQL(1)];
@testset "Testing Stochastic forcing structure with zero driving..."
        for eq in eqs
                sol = integrate(prob,eq,tspan;tsargs...);
                E = energy.(Ref(prob.d),sol.u)
                Z = enstrophy.(Ref(prob.d),sol.u)
                @test E[1] ≈ E[end] atol=1e-3
                @test Z[1] ≈ Z[end] atol=1e-3
        end
end
