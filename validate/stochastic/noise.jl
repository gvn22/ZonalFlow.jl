using Revise
using ZonalFlow

tspan   = (0.0,500.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=50,
            save_noise=true
           );

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=2,ε=0.005,τ=0.0,isotropic=false);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "validate/stochastic/8x8/constantinou/"

eqs     = [NL(),GQL(Λ=0),CE2(),GCE2(0),GQL(1),GCE2(1),GQL(2),GCE2(2)];
labels  = ["nl_gaussian","ql_gaussian","ce2_gaussian","gce2_0_gaussian","gql_1_gaussian","gce2_1_gaussian","gql_2_gaussian","gce2_2_gaussian"];

nl      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],nl,dn=dn,fn=labels[1])

ql      = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ql,dn=dn,fn=labels[2])

ce2      = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2,dn=dn,fn=labels[3])

gql     = integrate(prob,eqs[5],tspan;tsargs...)
write(prob,eqs[5],gql,dn=dn,fn=labels[5])

gce2      = integrate(prob,eqs[6],tspan;u0=gql.u[end],tsargs...)
write(prob,eqs[6],gce2,dn=dn,fn=labels[6])

gql_2     = integrate(prob,eqs[7],tspan;tsargs...)
write(prob,eqs[7],gql_2,dn=dn,fn=labels[7])

gce2_2      = integrate(prob,eqs[8],tspan;u0=gql_2.u[end],tsargs...)
write(prob,eqs[8],gce2_2,dn=dn,fn=labels[8])
