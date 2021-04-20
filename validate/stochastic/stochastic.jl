using Revise
using ZonalFlow

domain  = Domain(extent=(2.0π,2.0π),res=(6,6));
coeffs  = Coefficients(Ω=2.0π,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=true);
forcing = Stochastic(kf=3,dk=1,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,10000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=20
           );

eqs     = [GQL(Λ=0),CE2()];
labels  = ["ql_linear","ce2_linear"];
sols    = integrate(prob,eqs,tspan;tsargs...);
write(prob,eqs,sols,dn="validate/stochastic/6x6/",labels=labels)

eqs     = [NL(),GCE2(Λ=prob.d.nx-1)];
labels  = ["nl_linear","gce2_linear"];
write(prob,eqs,sols,dn="validate/stochastic/6x6/",labels=labels)
