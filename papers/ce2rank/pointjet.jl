using Revise
using ZonalFlow
using JLD2

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=10000,
            save_everystep=false,
            save_start=true,
            dense=false,
            save_noise=false,
            saveat=5
           );

eqs = [GQL(0),CE2(),CE2(),CE2()];
labels = ["ql_pj","ce2_pj","ce2_qlic_pj","ce2_fr_pj"];

dn      = "data/8x8/tau20/"
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

ce2_fr     = integrate(prob,eqs[4],tspan;u0=1e-9,tsargs...)
write(prob,eqs[4],ce2_fr,dn=dn,fn=labels[4])

nl      = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),nl,dn=dn,fn="nl")

dn      = "data/12x12/tau20/"
domain  = Domain(extent=(2π,2π),res=(12,12));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic= integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

ce2_fr  = integrate(prob,eqs[4],tspan;u0=1e-9,tsargs...)
write(prob,eqs[4],ce2_fr,dn=dn,fn=labels[4])

nl      = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),nl,dn=dn,fn="nl")

dn      = "data/16x16/tau20/"
domain  = Domain(extent=(2π,2π),res=(16,16));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic= integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

ce2_fr  = integrate(prob,eqs[4],tspan;u0=1e-9,tsargs...)
write(prob,eqs[4],ce2_fr,dn=dn,fn=labels[4])

nl      = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),nl,dn=dn,fn="nl")

dn      = "data/16x16/beta0/"
domain  = Domain(extent=(2π,2π),res=(16,16));
coeffs  = Coefficients(Ω=2π,θ=90.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic= integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...)
write(prob,eqs[3],ce2_qlic,dn=dn,fn="ce2_qlic_pj")

ce2_fr  = integrate(prob,eqs[4],tspan;u0=1e-9,tsargs...)
write(prob,eqs[4],ce2_fr,dn=dn,fn=labels[4])

nl      = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),nl,dn=dn,fn="nl")
