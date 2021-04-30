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

domain  = Domain(extent=(2π,2π),res=(12,12));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

eqs     = [NL(),GQL(Λ=prob.d.nx-1),GCE2(Λ=prob.d.nx-1),GQL(0),CE2()];
labels  = ["nl","gql","gce2","ql","ce2"];
sols    = integrate(prob,eqs,tspan;tsargs...);

dn      = "validate/pointjet/12x12/"

write(prob,eqs,ql,dn=dn,labels=labels)
