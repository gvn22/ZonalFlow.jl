using Revise
using ZonalFlow
using JLD2

tspan   = (0.0,5000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=50,
            save_noise=false
           );

# 6x6
# Couldn't find jets <- res possibly too small

# 8x8
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=2,ε=0.005);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "validate/stochastic/8x8/noisetests/"

eqs     = [NL(),GQL(Λ=0),CE2(),GCE2(0)];
labels  = ["nl","ql","ce2","gce2_0"];

nl      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],nl,dn=dn,fn=labels[1])

ql      = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ql,dn=dn,fn=labels[2])

ce2      = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2,dn=dn,fn=labels[3])

gce2      = integrate(prob,eqs[4],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[4],gce2,dn=dn,fn=labels[4])

# sols    = integrate(prob,eqs,tspan;tsargs...)
# write(prob,eqs,sols,dn=dn,labels=labels)

eqs     = [GQL(Λ=1),GCE2(Λ=1)];
labels  = ["gql_1_m2","gce2_1_m2"];

gql     = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],gql,dn=dn,fn=labels[1])

gce2    = integrate(prob,eqs[2],tspan;u0=gql.u[end],tsargs...)
write(prob,eqs[2],gce2,dn=dn,fn=labels[2])
@save dn*labels[2]*".jld2" gce2

eqs     = [GQL(Λ=2),GCE2(Λ=2)];
labels  = ["gql_2_m2","gce2_2_m2"];

gql     = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],gql,dn=dn,fn=labels[1])

tspan   = (0.0,500.0);
gce2    = integrate(prob,eqs[2],tspan;u0=gql.u[end],tsargs...)
write(prob,eqs[2],gce2,dn=dn,fn=labels[2])

eqs     = [GQL(Λ=3),GCE2(Λ=3)];
labels  = ["gql_3_m2","gce2_3_m2"];

tspan   = (0.0,2000.0);
gql     = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],gql,dn=dn,fn=labels[1])

tspan   = (0.0,500.0);
gce2    = integrate(prob,eqs[2],tspan;u0=gql.u[end],tsargs...)
write(prob,eqs[2],gce2,dn=dn,fn=labels[2])

# sols    = integrate(prob,eqs,tspan;tsargs...)
# write(prob,eqs,sols,dn=dn,labels=labels)

# 12x12
domain  = Domain(extent=(2π,2π),res=(12,12));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=8,dk=2,ε=0.001);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "validate/stochastic/12x12/"

eqs     = [GQL(Λ=0),CE2()];
labels  = ["ql_long","ce2_long"];

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])
@save dn*labels[1]*".jld2" ql

ce2     = integrate(prob,CE2(),tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])
@save dn*labels[2]*".jld2" ce2
