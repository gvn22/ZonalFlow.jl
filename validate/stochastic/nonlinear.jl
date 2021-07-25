using Revise
using ZonalFlow

tspan = (0.0,1000.0);
tsargs  = ( dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            save_noise=false,
            saveat=10
           );

eqs     = [GQL(Λ=0),CE2(),CE2()];

"""
   Single mode excitation
"""
dn      = "validate/stochastic/nonlinear/1.0/"

coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=3,dk=0,ε=1.0); #note dk = 0 => only kf is excited

labels  = ["ql_6x6_m3","ce2_6x6_m3","ce2_6x6_qlic_m3"];

domain  = Domain(extent=(2π,2π),res=(6,6));
prob    = BetaPlane(domain,coeffs,forcing);

ql     = integrate(prob,eqs[1],(0.0,5000.0);tsargs...);
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...);
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

labels  = ["ql_8x8_m3","ce2_8x8_m3","ce2_8x8_qlic_m3"];

domain  = Domain(extent=(2π,2π),res=(8,8));
prob    = BetaPlane(domain,coeffs,forcing);

ql     = integrate(prob,eqs[1],(0.0,5000.0);tsargs...);
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...);
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

labels  = ["ql_10x10_m3","ce2_10x10_m3","ce2_10x10_qlic_m3"];

domain  = Domain(extent=(2π,2π),res=(10,10));
prob    = BetaPlane(domain,coeffs,forcing);

ql     = integrate(prob,eqs[1],(0.0,5000.0);tsargs...);
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...);
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

"""
   Three-mode excitations
"""
dn      = "validate/stochastic/nonlinear/e0.1m5d1/"

coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=1,ε=0.1);

# 6x6

labels  = ["ql_6x6","ce2_6x6","ce2_6x6_qlic"];

domain  = Domain(extent=(2π,2π),res=(6,6));
prob    = BetaPlane(domain,coeffs,forcing);

ql     = integrate(prob,eqs[1],(0.0,5000.0);tsargs...);
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...);
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

# 8x8

domain  = Domain(extent=(2π,2π),res=(8,8));
prob    = BetaPlane(domain,coeffs,forcing);

labels  = ["ql_8x8","ce2_8x8","ce2_8x8_qlic"];

ql     = integrate(prob,eqs[1],(0.0,10000.0);tsargs...);
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...);
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

# 10x10

domain  = Domain(extent=(2π,2π),res=(10,10));
prob    = BetaPlane(domain,coeffs,forcing);

labels  = ["ql_10x10","ce2_10x10","ce2_10x10_qlic"];

ql     = integrate(prob,eqs[1],(0.0,5000.0);tsargs...);
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,CE2(poscheckat=10),tspan;u0=ql.u[end],tsargs...);
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

# 12x12

domain  = Domain(extent=(2π,2π),res=(12,12));
prob    = BetaPlane(domain,coeffs,forcing);

labels  = ["ql_12x12","ce2_12x12","ce2_12x12_qlic"];

ql     = integrate(prob,eqs[1],(0.0,5000.0);tsargs...);
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

ce2_qlic     = integrate(prob,eqs[3],tspan;u0=ql.u[end],tsargs...);
write(prob,eqs[3],ce2_qlic,dn=dn,fn=labels[3])

# GCE2
# eqs     = [NL(),GQL(Λ=prob.d.nx-1),GCE2(Λ=prob.d.nx-1)];
# labels  = ["nl_linear","gql_5_linear","gce2_5_linear"];
# sols    = integrate(prob,eqs,tspan;tsargs...);
# write(prob,eqs,sols,dn="validate/stochastic/6x6/",labels=labels)
