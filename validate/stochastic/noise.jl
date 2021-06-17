using Revise
using ZonalFlow
using JLD2

tspan   = (0.0,20000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=50000,
            save_everystep=false,
            saveat=50,
            dense=false,
            save_noise=false
           );

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=2,ε=0.005,τ=0.0,isotropic=false);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "validate/stochastic/8x8/bandtype/"

eqs     = [NL(),GQL(7),GCE2(7),GQL(Λ=0),CE2(),GCE2(0),GQL(1),GCE2(1),GQL(2),GCE2(2)];
labels  = ["nl","gql_m","gce2_m","ql","ce2","gce2_0","gql_1","gce2_1","gql_2","gce2_2"];

nl      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],nl,dn=dn,fn=labels[1])
@save dn*labels[1]*".jld2" nl

gql_7      = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],gql_7,dn=dn,fn=labels[2])
@save dn*labels[2]*".jld2" gql_7

gce2_7     = integrate(prob,eqs[3],tspan;tsargs...)
write(prob,eqs[3],gce2_7,dn=dn,fn=labels[3])
@save dn*labels[3]*".jld2" gce2_7

ql      = integrate(prob,eqs[4],tspan;tsargs...)
write(prob,eqs[4],ql,dn=dn,fn="ql_zero")
@save dn*labels[4]*".jld2" ql

ce2      = integrate(prob,eqs[5],tspan;tsargs...)
write(prob,eqs[5],ce2,dn=dn,fn=labels[5])
@save dn*labels[5]*".jld2" ce2

gce2      = integrate(prob,eqs[6],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[6],gce2,dn=dn,fn=labels[6])
@save dn*labels[6]*".jld2" gce2

gql_1     = integrate(prob,eqs[7],tspan;tsargs...)
write(prob,eqs[7],gql_1,dn=dn,fn=labels[7])
@save dn*labels[7]*".jld2" gql_1

gce2_1    = integrate(prob,eqs[8],tspan;u0=gql_1.u[end],tsargs...)
write(prob,eqs[8],gce2_1,dn=dn,fn=labels[8])
@save dn*labels[8]*".jld2" gce2_1
# @load dn*labels[8]*".jld2" gce2_1

gql_2     = integrate(prob,eqs[9],tspan;tsargs...)
write(prob,eqs[9],gql_2,dn=dn,fn="gql_2")
@save dn*labels[9]*".jld2" gql_2

tspan = (0.0,2000.0)
gce2_2     = integrate(prob,GCE2(2),tspan;tsargs...)
write(prob,GCE2(2),gce2_2,dn=dn,fn="gce2_2_zeros")
@save dn*labels[10]*"_zeros.jld2" gce2_2
