using Revise
using ZonalFlow

domain  = Domain(extent=(2π,2π),res=(6,6));
tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=20
           );

# Nonlinear
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=3,dk=1,ε=0.001);
prob    = BetaPlane(domain,coeffs,forcing);

eqs     = [GQL(Λ=0),CE2()];
labels  = ["ql","ce2"];
sols    = integrate(prob,eqs,tspan;tsargs...);
write(prob,eqs,sols,dn="validate/stochastic/6x6/",labels=labels)

eqs     = [NL(),GQL(Λ=prob.d.nx-1),GCE2(Λ=prob.d.nx-1)];
labels  = ["nl","gql_5","gce2_5"];
sols    = integrate(prob,eqs,tspan;tsargs...);
write(prob,eqs,sols,dn="validate/stochastic/6x6/",labels=labels)

# Linear
domain  = Domain(extent=(2π,2π),res=(6,6));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=true);
forcing = Stochastic(kf=3,dk=1,ε=0.001);
prob    = BetaPlane(domain,coeffs,forcing);

eqs     = [GQL(Λ=0),CE2()];
labels  = ["ql_linear","ce2_linear"];
sols    = integrate(prob,eqs,tspan;tsargs...);
write(prob,eqs,sols,dn="validate/stochastic/6x6/",labels=labels)

eqs     = [NL(),GQL(Λ=prob.d.nx-1),GCE2(Λ=prob.d.nx-1)];
labels  = ["nl_linear","gql_5_linear","gce2_5_linear"];
sols    = integrate(prob,eqs,tspan;tsargs...);
write(prob,eqs,sols,dn="validate/stochastic/6x6/",labels=labels)
