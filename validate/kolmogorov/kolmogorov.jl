using Revise
using ZonalFlow
using JLD2

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=5
           );

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

eqs     = [NL(),GQL(7),GCE2(7),GQL(0),CE2(),GCE2(Λ=0,poscheck=true,poscheckat=10),GQL(1),GCE2(1),GQL(3),GCE2(3)];
labels  = ["nl","gql_m","gce2_m","ql_ts5","ce2_ts5","gce2_0","gql_1","gce2_1","gql_3","gce2_3"];
# sols    = integrate(prob,eqs,tspan;tsargs...);

dn      = "validate/kolmogorov/8x8/"

nl      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],nl,dn=dn,fn=labels[1])
@save dn*labels[1]*".jld2" nl

gql_1      = integrate(prob,GQL(1),tspan;tsargs...)

gql_7      = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],gql_7,dn=dn,fn=labels[2])
@save dn*labels[2]*".jld2" gql_7

gce2_7     = integrate(prob,eqs[3],tspan;tsargs...)
write(prob,eqs[3],gce2_7,dn=dn,fn=labels[3])
@save dn*labels[3]*".jld2" gce2_7

ql      = integrate(prob,eqs[4],tspan;u0=0.0,tsargs...)
write(prob,eqs[4],ql,dn=dn,fn=labels[4])
@save dn*labels[4]*".jld2" ql
@load dn*labels[4]*".jld2" ql

ce2     = integrate(prob,eqs[5],tspan;tsargs...)
write(prob,eqs[5],ce2,dn=dn,fn=labels[5])
@save dn*labels[5]*".jld2" ce2

gce2      = integrate(prob,eqs[6],tspan;tsargs...)
write(prob,eqs[6],gce2,dn=dn,fn=labels[6])
@save dn*labels[6]*".jld2" gce2

gql_1     = integrate(prob,eqs[7],tspan;tsargs...)
write(prob,eqs[7],gql_1,dn=dn,fn=labels[7])
@save dn*labels[7]*".jld2" gql_1

gce2_1    = integrate(prob,eqs[8],tspan;u0=gql_1.u[end],tsargs...)
write(prob,eqs[8],gce2_1,dn=dn,fn=labels[8])
@save dn*labels[8]*".jld2" gce2_1
# @load dn*labels[8]*".jld2" gce2_1

@info "Running:" eqs[9]
gql_3     = integrate(prob,eqs[9],tspan;tsargs...)
write(prob,eqs[9],gql_2,dn=dn,fn=labels[9])
@save dn*labels[9]*".jld2" gql_2

@info "Running:" eqs[10]
gce2_3     = integrate(prob,eqs[10],tspan;tsargs...)
write(prob,eqs[10],gce2_3,dn=dn,fn=labels[10])
@save dn*labels[10]*".jld2" gce2_3

# sol    = integrate(prob,NL(),tspan;tsargs...);
# write(prob,NL(),sol,dn="validate/kolmogorov/6x6/",fn="nl")
#
# sol    = integrate(prob,GQL(Λ=prob.d.nx-1),tspan;tsargs...);
# write(prob,GQL(Λ=prob.d.nx-1),sol,dn="validate/kolmogorov/6x6/",fn="gql")
#
# sol    = integrate(prob,GCE2(Λ=prob.d.nx-1),tspan;tsargs...);
# write(prob,GCE2(Λ=prob.d.nx-1),sol,dn="validate/kolmogorov/6x6/",fn="gce2")
#
# sol    = integrate(prob,GQL(0),tspan;tsargs...);
# write(prob,GQL(0),sol,dn="validate/kolmogorov/6x6/",fn="ql")
#
# sol    = integrate(prob,CE2(),tspan;u0=sol.u[end],tsargs...);
# write(prob,CE2(),sol,dn="validate/kolmogorov/6x6/",fn="ce2")
#
# sol    = integrate(prob,GCE2(0),tspan;tsargs...);
# write(prob,GCE2(0),sol,dn="validate/kolmogorov/6x6/",fn="gce2_0")
#
# sol    = integrate(prob,GQL(1),tspan;tsargs...);
# write(prob,GQL(0),sol,dn="validate/kolmogorov/6x6/",fn="gql_1")
#
# sol    = integrate(prob,GCE2(1),tspan;tsargs...);
# write(prob,GQL(0),sol,dn="validate/kolmogorov/6x6/",fn="gce2_1")
