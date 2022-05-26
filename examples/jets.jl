using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using ZonalFlow

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=100,
            save_noise=false
           );

domain  = Domain(extent=(2π,2π),res=(16,16));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0);
forcing = Stochastic(kf=10,dk=4,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);
eq      = NL()

sol = integrate(prob,eq,tspan;tsargs...);
write(prob,eq,sol,dn="data/",fn="jets_nl")
