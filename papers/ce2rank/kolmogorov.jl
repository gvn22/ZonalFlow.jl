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
            save_start=true,
            saveat=5
           );

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "data/8x8/rich/"

eqs = [GQL(0),CE2(),CE2(),CE2(eigmax=true)];
labels = ["ql_kf_zm","ce2_kf","ce2_qlic_kf","ce2_qlic_eigmax_kf"];

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])
@save "ql_kf_.jld2" ql

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

ce2_qlic_eigmax     = integrate(prob,eqs[4],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[4],ce2_qlic_eigmax,dn=dn,fn=labels[4])

ce2_qlic_m1     = integrate(prob,CE2(),tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic_m1,dn=dn,fn="ce2_qlic_m1r3t10_av_kf")

ce2_qlic_m2     = integrate(prob,CE2(),tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic_m2,dn=dn,fn="ce2_qlic_m2_kf")

ql_qlic      = integrate(prob,eqs[1],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[1],ql_qlic,dn=dn,fn="ql_qlic_kf")

@load dn*"ql.jld2" ql
ql_qlic_m2      = integrate(prob,eqs[1],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[1],ql_qlic_m2,dn=dn,fn="ql_qlic_m2_kf_2")

# dn      = "data/8x8/brad/"

# nl      = integrate(prob,NL(),tspan;tsargs...)
# write(prob,NL(),nl,dn=dn,fn="nl_kf")

domain  = Domain(extent=(2π,2π),res=(8,12));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "data/8x12/nu02/"

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

ce2_qlic_eigmax     = integrate(prob,eqs[4],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[4],ce2_qlic_eigmax,dn=dn,fn=labels[4])

domain  = Domain(extent=(2π,2π),res=(12,12));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "data/12x12/nu02/"

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

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,CE2(poscheck=true,poscheckat=1),tspan;u0=1e-16,tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn="ce2_kf_2")

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

ce2_qlic_eigmax     = integrate(prob,eqs[4],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[4],ce2_qlic_eigmax,dn=dn,fn=labels[4])
