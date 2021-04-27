using Revise
using ZonalFlow

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=50
           );

domain  = Domain(extent=(2π,2π),res=(6,6));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

# eqs     = [NL(),GQL(Λ=prob.d.nx-1),GCE2(Λ=prob.d.nx-1),QL(),CE2(),GCE2(Λ=0)];
# labels  = ["nl","gql","gce2","ql","ce2","gce2_0"];
# sols    = integrate(prob,eqs,tspan;tsargs...);

sol    = integrate(prob,NL(),tspan;tsargs...);
write(prob,NL(),sol,dn="validate/kolmogorov/6x6/",fn="nl")

sol    = integrate(prob,GQL(Λ=prob.d.nx-1),tspan;tsargs...);
write(prob,GQL(Λ=prob.d.nx-1),sol,dn="validate/kolmogorov/6x6/",fn="gql")

sol    = integrate(prob,GCE2(Λ=prob.d.nx-1),tspan;tsargs...);
write(prob,GCE2(Λ=prob.d.nx-1),sol,dn="validate/kolmogorov/6x6/",fn="gce2")

sol    = integrate(prob,GQL(0),tspan;tsargs...);
write(prob,GQL(0),sol,dn="validate/kolmogorov/6x6/",fn="ql")

sol    = integrate(prob,CE2(),tspan;u0=sol.u[end],tsargs...);
write(prob,CE2(),sol,dn="validate/kolmogorov/6x6/",fn="ce2")

sol    = integrate(prob,GCE2(0),tspan;tsargs...);
write(prob,GCE2(0),sol,dn="validate/kolmogorov/6x6/",fn="gce2_0")

sol    = integrate(prob,GQL(1),tspan;tsargs...);
write(prob,GQL(0),sol,dn="validate/kolmogorov/6x6/",fn="gql_1")

sol    = integrate(prob,GCE2(1),tspan;tsargs...);
write(prob,GQL(0),sol,dn="validate/kolmogorov/6x6/",fn="gce2_1")
