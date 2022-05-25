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

f1 = PointJet(Ξ=0.0,Δθ=0.1,τ=10.0);
f2 = Kolmogorov(A₁=0.0,A₄=0.0);
f3 = Stochastic(kf=3,dk=1,ε=0.0);
fls = ["Point Jet", "Kolmogorov", "Stochastic"]

eqs = [NL(),CE2(),GQL(0)];
# eqs     = append!(eqs,[GQL(Λ=l) for l=0:prob.d.nx-1])
# eqs     = append!(eqs,[GCE2(Λ=l) for l=0:prob.d.nx-1])

for (label,forcing) in zip(fls,[f1,f2,f3])
        prob = BetaPlane(domain,coeffs,forcing);
        for eq in eqs
                sol = integrate(prob,eq,tspan;tsargs...);
                E = energy.(prob.d,sol.u)
                Z = enstrophy.(prob.d,sol.u)
                @test E[1] == E[end]
                @test Z[1] == Z[end]
        end
end
