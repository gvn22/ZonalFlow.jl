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
            saveat=10
           );

eqs = [GQL(0),CE2(),CE2()];
labels = ["ql","ce2","ce2_qlic"];

dn      = "data/8x8/sftests_sriw_conservative/"
# test conservative
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0,linear=true);
forcing = Stochastic(kf=5,dk=0,ε=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,1000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],(0.0,1000.0);tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_nl_conservative/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,1000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],(0.0,1000.0);tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_nl_o1_conservative/"

tsargs  = (
            dt=0.0005,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=10
           );

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1]*"dt0005")

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sosra_conservative/"
# test conservative
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0,linear=true);
forcing = Stochastic(kf=5,dk=0,ε=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sosra_linear/"
# test conservative
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=true);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,4000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1]*"nu4")

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=10
           );

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2]*"nu4")

dn      = "data/8x8/sftests_sosra_nonlinear/"
# test conservative
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=10
           );

ql      = integrate(prob,eqs[1],(0.0,4000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1]*"nu4")

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2]*"nu4")

dn      = "data/8x8/sftests_sriw1_nonlinear/"
# test conservative
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=10
           );

ql      = integrate(prob,eqs[1],(0.0,4000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1]*"nu4")

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2]*"nu4")

dn      = "data/8x8/sftests_sriw1_nonlinear_0m3/"
# test conservative
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=10
           );

ql      = integrate(prob,eqs[1],(0.0,4000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1]*"nu4")

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2]*"nu4")

dn      = "data/8x8/sftests_sriw_nl_conservative_m0/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_m0/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_noise/"

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=true,
            saveat=10
           );

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,200.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],(0.0,2500.0);tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_m0_r1e-4/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,20000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],(0.0,2500.0);tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_m0_r1e-4_m30_dt001/tests/"

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=10
           );

tspan = (0.0,20000.0);
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql_m30_dt0005      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql_m30_dt0005,dn=dn,fn=labels[1])
ql = ql_m30_dt0005
@save dn*"ql.jld2" ql

ce2     = integrate(prob,eqs[2],(0.0,5000.0);tsargs...)
@save dn*"ce2.jld2" ce2
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sosra_m0_r1e-4_m30_dt0005/"

tsargs  = (
            dt=0.0005,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=10
           );

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql_m30_dt0005_sosra      = integrate(prob,eqs[1],(0.0,25000.0);tsargs...)
write(prob,eqs[1],ql_m30_dt0005_sosra,dn=dn,fn=labels[1]*"test")

u0 = deepcopy(ql_m30_dt0005_sosra.u[end])

for m = 1:7
    for n = -7:7
        u0[n+prob.d.ny,m+1] = 1e-4*exp(im*rand(Uniform(0,2π)))
    end
end

u0[:,1]

ce2_m30     = integrate(prob,eqs[2],(0.0,1000.0);u0=u0,tsargs...)
write(prob,eqs[2],ce2_m30,dn=dn,fn=labels[2]*"dt001")

@info length(ql_m30_dt0005_sosra.u)
u0 = deepcopy(mean(ql_m30_dt0005_sosra.u[2000]))
u0[:,1]
u0 = ql_m30_dt0005_sosra.u[2000]
ce2_m30_av     = integrate(prob,eqs[2],(0.0,1000.0);u0=u0,tsargs...)
write(prob,eqs[2],ce2_m30_av,dn=dn,fn=labels[2]*"dt001av")

dn      = "data/8x8/sftests_sriw_m0_r1e-4_m30/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql_m30      = integrate(prob,eqs[1],(0.0,25000.0);tsargs...)
write(prob,eqs[1],ql_m30,dn=dn,fn=labels[1])

ce2_m30     = integrate(prob,eqs[2],(0.0,2500.0);tsargs...)
write(prob,eqs[2],ce2_m30,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_m0_r1e-5/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,10000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_m0_r1e-6/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,10000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_qlce2_fp/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,10000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_m0_r1e-3/"

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            dense=false,
            save_noise=false,
            saveat=10
           );

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,1000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])


dn      = "data/8x8/sftests_sriw_m0_mu/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=1.0,ν=0.001,ν₄=1.0,linear=true);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_m0_mu2/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=100.0,ν=0.001,ν₄=1.0,linear=true);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/sftests_sriw_m0_m30/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.01,ν=0.001,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=0,ε=0.05);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],(0.0,4000.0);tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/pjtests_conservative/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0,linear=true);
forcing = PointJet(Ξ=0.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

eqs = [GQL(0),CE2(),CE2()];
labels = ["ql","ce2","ce2_qlic"];

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/kftests_conservative/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0,linear=true);
forcing = Kolmogorov(A₁=0.0,A₄=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/kftests_conservative_rand0_1/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=8.0,θ=0.0,μ=0.0,ν=0.0,ν₄=0.0,linear=true);
forcing = Kolmogorov(A₁=0.0,A₄=0.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/kftests_nonlinear_rand0_1/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/pjtests_m0/"
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])

dn      = "data/8x8/kftests_m0/"

domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql,dn=dn,fn=labels[1])

ce2     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2,dn=dn,fn=labels[2])
