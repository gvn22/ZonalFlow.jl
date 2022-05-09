using Revise
using ZonalFlow
using JLD2

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=200000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=50
           );

eqs = [GQL(0),CE2(),CE2(),GQL(1),GCE2(1)];
labels = ["ql","ce2","ce2_qlic","gql_1","gce2_1"];

dn      = "data/8x8/nu01/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.005);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,2000.0);

nl      = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),nl,dn=dn,fn="nl")

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

dn      = "data/8x8/nu01e01/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,2000.0);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

tspan   = (0.0,5000.0);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1]*"t5k")

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2]*"t5k")

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3]*"t5k")

ce2_fr     = integrate(prob,eqs[3],tspan;u0=1e-6,tsargs...)
write(prob,eqs[3],ce2_fr,dn=dn,fn="ce2frt5k")

dn      = "data/8x8/mu01e01m04n023/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=4,dk=0,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,1000.0);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

dn      = "data/8x8/nu01e005k04/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=4,dk=0,ε=0.005);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,2000.0);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

dn      = "data/8x8/mu01e01m04n03_2/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=4,dk=0,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,1000.0);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

dn      = "data/8x8/mu01e01m04n02/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=4,dk=0,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,1000.0);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])
