using Revise
using Traceur
using BenchmarkTools
using ZonalFlow

domain  = Domain(extent=(2π,2π),res=(4,4));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.01,ν=1.0,ν₄=1.0,linear=false);
forcing = Stochastic(kf=3,dk=1,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            save_start=false
           );

# eq = NL()
# eq = CE2()
# eq = GQL(2)
eq = GCE2(2)

u0 = get_de_ic(prob,eq);
p  = get_de_params(prob,eq)

du = similar(u0);
@benchmark f!($du,$u0,$p,$tspan)

_prob,_alg = get_de_probalg(prob,eq,u0,tspan,p);
@benchmark solve($_prob,$_alg;$tsargs...)
