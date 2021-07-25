using Revise
using ZonalFlow
using JLD2

tspan   = (0.0,50.0);
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

dn      = "data/8x8/"

eqs = [GQL(0),CE2(eigmax=true)];
labels = ["ql_test","ce2_kf_eigmax_qlic_cb"];

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])
