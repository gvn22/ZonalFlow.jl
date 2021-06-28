using Revise
using ZonalFlow
using JLD2

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=50000,
            save_everystep=false,
            save_start=true,
            # dense=true,
            save_noise=false,
            saveat=5
           );

dn      = "validate/ranktests/repeat/";
eqs     = [GQL(0),CE2()];

## Point Jet
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql_pj      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql_pj,dn=dn,fn="ql_pj_500")

ce2_pj     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2_pj,dn=dn,fn="ce2_pj_500")

ce2_pj     = integrate(prob,eqs[2],tspan;u0=ql_pj.u[end],tsargs...)
write(prob,eqs[2],ce2_pj,dn=dn,fn="ce2_pj_qlic_500")

## Kolmogorov forcing
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql_kf      = integrate(prob,eqs[1],tspan;tsargs...);
write(prob,eqs[1],ql_kf,dn=dn,fn="ql_kf_500");

ce2_kf     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2_kf,dn=dn,fn="ce2_kf_500");

ce2_kf_qlic     = integrate(prob,eqs[2],tspan;u0=ql_kf.u[end],tsargs...);
write(prob,eqs[2],ce2_kf_qlic,dn=dn,fn="ce2_kf_qlic_500");

"""
# Kolmogorov with positivity check
eqs     = [GQL(0),CE2(poscheck=true,poscheckat=10)];

ce2_kf     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2_kf,dn=dn,fn="ce2_kf_poscheck_500");

ce2_kf_qlic     = integrate(prob,eqs[2],tspan;u0=ql_kf.u[end],tsargs...);
write(prob,eqs[2],ce2_kf_qlic,dn=dn,fn="ce2_kf_poscheck_qlic_500");
"""

"""
# # Kolmogorov forcing lower viscosity
domain  = Domain(extent=(2π,2π),res=(11,11));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.01,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

eqs     = [GQL(0),CE2(poscheck=true,poscheckat=5)];

ql_kf      = integrate(prob,eqs[1],tspan;tsargs...);
write(prob,eqs[1],ql_kf,dn=dn,fn="ql_kf_nu01_500");

ce2_kf     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2_kf,dn=dn,fn="ce2_kf_nu01_500");

ce2_kf_qlic     = integrate(prob,eqs[2],tspan;u0=ql_kf.u[end],tsargs...);
write(prob,eqs[2],ce2_kf_qlic,dn=dn,fn="ce2_kf_nu01_qlic_500");
"""

# Kolmogorov Forcing: 11x11 grid
domain  = Domain(extent=(2π,2π),res=(9,9));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

ql_kf_9      = integrate(prob,eqs[1],tspan;tsargs...);
write(prob,eqs[1],ql_kf_9,dn=dn,fn="ql_kf_9x9_500");

ce2_kf_9     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2_kf_9,dn=dn,fn="ce2_kf_9x9_500");

ce2_kf_9_qlic     = integrate(prob,eqs[2],tspan;u0=ql_kf_9.u[end],tsargs...);
write(prob,eqs[2],ce2_kf_9_qlic,dn=dn,fn="ce2_kf_9x9_qlic_500");

# Kolmogorov Forcing: 11x11 grid
domain  = Domain(extent=(2π,2π),res=(11,11));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,1000.0);

ql_kf_11      = integrate(prob,eqs[1],tspan;tsargs...);
write(prob,eqs[1],ql_kf_11,dn=dn,fn="ql_kf_11x11_1000");

tsargs  = (
            dt=0.0005,
            adaptive=false,
            progress=true,
            progress_steps=50000,
            save_everystep=false,
            save_start=true,
            # dense=true,
            save_noise=false,
            saveat=5
           );

ce2_kf_11     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2_kf_11,dn=dn,fn="ce2_kf_11x11_500");
@save dn*"ce2_kf_11x11_500"*".jld2" ce2_kf_11

ce2_kf_11_qlic     = integrate(prob,eqs[2],tspan;u0=ql_kf_11.u[end],tsargs...);
write(prob,eqs[2],ce2_kf_11_qlic,dn=dn,fn="ce2_kf_11x11_qlic_500");
@save dn*"ce2_kf_qlic_11x11_500"*".jld2" ce2_kf_11_qlic

# Kolmogorov Forcing: 16x16 grid
domain  = Domain(extent=(2π,2π),res=(16,16));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

eqs     = [GQL(0),CE2(poscheck=true,poscheckat=5)];

tsargs  = (
            dt=0.0005,
            adaptive=false,
            progress=true,
            progress_steps=50000,
            save_everystep=false,
            save_start=true,
            # dense=true,
            save_noise=false,
            saveat=5
           );

ql_kf_16      = integrate(prob,eqs[1],tspan;tsargs...);
write(prob,eqs[1],ql_kf_16,dn=dn,fn="ql_kf_16x16_1000");
@save dn*"ql_kf_16x16_500"*".jld2" ql_kf_16

ce2_kf_16     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2_kf_16,dn=dn,fn="ce2_kf_16x16_500");
@save dn*"ce2_kf_16x16_500"*".jld2" ce2_kf_16

ce2_kf_16_qlic     = integrate(prob,eqs[2],tspan;u0=ql_kf_16.u[end],tsargs...);
write(prob,eqs[2],ce2_kf_16_qlic,dn=dn,fn="ce2_kf_16x16_qlic_500");
@save dn*"ce2_kf_16x16_qlic_500"*".jld2" ce2_kf_16_qlic

# Kolmogorov Forcing: 24x24 grid
domain  = Domain(extent=(2π,2π),res=(24,24));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=0.0); # β=0
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);

eqs     = [GQL(0),CE2()];

tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=5000,
            save_everystep=false,
            save_start=true,
            # dense=true,
            save_noise=false,
            saveat=20
           );

ql_kf_24      = integrate(prob,eqs[1],tspan;tsargs...);
write(prob,eqs[1],ql_kf_24,dn=dn,fn="ql_kf_24x24_500");
@save dn*ql_kf_24x24_500*".jld2" ql_kf_24

ce2_kf_24     = integrate(prob,eqs[2],tspan;tsargs...);
write(prob,eqs[2],ce2_kf_24,dn=dn,fn="ce2_kf_24x24_500");

## Stochastic jet
domain  = Domain(extent=(2π,2π),res=(8,8));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=true);
forcing = Stochastic(kf=5,dk=2,ε=0.005,τ=0.0,isotropic=false);
prob    = BetaPlane(domain,coeffs,forcing);

ql_sf      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql_sf,dn=dn,fn="ql_sf_test_500")

ce2_sf     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2_sf,dn=dn,fn="ce2_sf_test_500")

ce2_sf_qlic     = integrate(prob,eqs[2],tspan;u0=ql_sf.u[end],tsargs...)
write(prob,eqs[2],ce2_sf_qlic,dn=dn,fn="ce2_sf_qlic_500")

#
domain  = Domain(extent=(2π,2π),res=(12,12));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=5,dk=2,ε=0.005,τ=0.0,isotropic=false);
prob    = BetaPlane(domain,coeffs,forcing);

ql_sf      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql_sf,dn=dn,fn="ql_sf_test_1000")

ce2_sf     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2_sf,dn=dn,fn="ce2_sf_test_1000")

tspan = (0.0,500.0)
domain  = Domain(extent=(2π,2π),res=(12,12));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=8,dk=2,ε=0.01,τ=0.0,isotropic=false);
prob    = BetaPlane(domain,coeffs,forcing);

ql_sf      = integrate(prob,eqs[1],tspan;tsargs...)
write(prob,eqs[1],ql_sf,dn=dn,fn="ql_sf_test_1000")

ce2_sf     = integrate(prob,eqs[2],tspan;tsargs...)
write(prob,eqs[2],ce2_sf,dn=dn,fn="ce2_sf_test_1000")

ce2_sf_qlic     = integrate(prob,eqs[2],tspan;u0=ql_sf.u[end],tsargs...)
write(prob,eqs[2],ce2_sf_qlic,dn=dn,fn="ce2_sf_qlic_test_1000")
