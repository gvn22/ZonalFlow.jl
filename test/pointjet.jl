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

domain  = Domain(extent=(2π,2π),res=(5,5));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=0.0,Δθ=0.1,τ=10.0);
prob    = BetaPlane(domain,coeffs,forcing);

eqs = [NL(),GQL(0),GQL(2)];
@testset "Testing Point Jet structure with zero driving..." begin
        for eq in eqs
                sol = integrate(prob,eq,tspan;tsargs...);
                E = energy.(Ref(prob.d),sol.u)
                Z = enstrophy.(Ref(prob.d),sol.u)
                @test E[1] ≈ E[end]
                @test Z[1] ≈ Z[end]
        end
end
