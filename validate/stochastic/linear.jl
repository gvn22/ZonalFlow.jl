using Revise
using ZonalFlow

domain  = Domain(extent=(2π,2π),res=(6,6));
tspan   = (0.0,10000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=50
           );

dn      = "validate/stochastic/6x6/"

domain  = Domain(extent=(2π,2π),res=(6,6));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=0.0,linear=true);
forcing = Stochastic(kf=3,dk=1,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

eqs     = [GQL(Λ=0),CE2()];
labels  = ["ql_linear","ce2_linear_qlic"];

sol     = integrate(prob,eqs[1],tspan;tsargs...);
write(prob,eqs[1],sol,dn=dn,fn=labels[1])

# use QL solution as u0
sol     = integrate(prob,eqs[2],tspan;u0=sol.u[end],tsargs...);
write(prob,eqs[2],sol,dn=dn,fn=labels[2])

eqs     = [NL(),GQL(Λ=prob.d.nx-1),GCE2(Λ=prob.d.nx-1)];
labels  = ["nl_linear","gql_5_linear","gce2_5_linear"];
sols    = integrate(prob,eqs,tspan;tsargs...);
write(prob,eqs,sols,dn="validate/stochastic/6x6/",labels=labels)
