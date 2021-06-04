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

domain  = Domain(extent=(2π,2π),res=(8,8));
# domain  = Domain(extent=(2π,2π),res=(12,12));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "validate/pointjet/8x8/"
# dn      = "validate/pointjet/12x12/"

eqs     = [NL(),GQL(Λ=prob.d.nx-1),GCE2(Λ=prob.d.nx-1),GQL(0),CE2()];
labels  = ["nl","gql","gce2","ql","ce2"];
sols    = integrate(prob,eqs,tspan;tsargs...);
write(prob,eqs,sols,dn=dn,labels=labels)

gce2_0    = integrate(prob,GCE2(0),tspan;tsargs...);
write(prob,GCE2(0),gce2_0,dn=dn,fn="gce2_0")

gql    = integrate(prob,GQL(1),tspan;tsargs...);
write(prob,GQL(1),gql,dn=dn,fn="gql_1")

gce2    = integrate(prob,GCE2(Λ=1,poscheck=true,poscheckat=20),tspan;tsargs...);
write(prob,GCE2(1),gce2,dn=dn,fn="gce2_1")

gql    = integrate(prob,GQL(3),tspan;tsargs...);
write(prob,GQL(3),gql,dn=dn,fn="gql_3")

gce2    = integrate(prob,GCE2(Λ=3,poscheck=true,poscheckat=20),tspan;tsargs...);
write(prob,GCE2(3),gce2,dn=dn,fn="gce2_3")

gql    = integrate(prob,GQL(5),tspan;tsargs...);
write(prob,GQL(5),gql,dn=dn,fn="gql_5")

gce2    = integrate(prob,GCE2(5),tspan;tsargs...);
write(prob,GCE2(5),gce2,dn=dn,fn="gce2_5")
