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
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=2,ε=0.005);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "data/8x8/stochastic/"

eqs = [GQL(0),CE2(),CE2(),CE2(eigmax=true)];
labels = ["ql_sf_test","ce2_sf_test","ce2_qlic_sf_test","ce2_qlic_eigmax_sf_test"];

nl      = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),nl,dn=dn,fn="nl_sf_call")

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=2,ε=0.005);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "data/8x8/stochastic/"

eqs = [GQL(0),CE2(),CE2(),CE2(eigmax=true)];
labels = ["ql_sf_c52","ce2_sf_c52","ce2_qlic_sf_c52","ce2_qlic_eigmax_sf_c52"];

nl      = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),nl,dn=dn,fn="nl_sf_c")

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

ce2_qlic_eigmax     = integrate(prob,eqs[4],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[4],ce2_qlic_eigmax,dn=dn,fn=labels[4])

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=4,dk=4,ε=0.005);
prob    = BetaPlane(domain,coeffs,forcing);

eqs = [GQL(0),CE2(),CE2(),CE2(eigmax=true)];
labels = ["ql_sf_c44","ce2_sf_c44","ce2_qlic_sf_c44","ce2_qlic_eigmax_sf_c44"];

nl      = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),nl,dn=dn,fn="nl_sf_c44")

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])
