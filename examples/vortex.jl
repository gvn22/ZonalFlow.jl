using ZonalFlow

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=10000,
            save_everystep=false,
            saveat=50,
            save_noise=false
           );

domain  = Domain(extent=(4π,2π),res=(24,24));
coeffs  = Coefficients(Ω=0.0,θ=0.0,μ=0.0,ν=0.02,ν₄=1.0);
forcing = Kolmogorov(A₁=-1.0,A₄=-2.0);
prob    = BetaPlane(domain,coeffs,forcing);
eq      = NL()

sol = integrate(prob,eq,tspan;tsargs...);
write(prob,eq,sol,dn="data/",fn="vortex_nl")
