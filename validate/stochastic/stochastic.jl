using Revise
using ZonalFlow

domain  = Domain(extent=(2.0π,2.0π),res=(6,6));
coeffs  = Coefficients(Ω=2.0π,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0,linear=true);
forcing = Stochastic(kf=3,dk=1,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

tspan   = (0.0,500.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=100000,
            save_everystep=false,
            saveat=100
           );

eqs     = [NL(),GQL(Λ=0),CE2()];
labels  = ["nl_linear","ql_linear","ce2_linear"];

@time sol1  = integrate(prob,eqs[1],tspan;tsargs...);
@time sol2  = integrate(prob,eqs[2],tspan;tsargs...);
@time sol3  = integrate(prob,eqs[3],tspan;tsargs...);
# @time sol4  = integrate(prob,GQL(Λ=1),tspan;tsargs...);
# @time sol5  = integrate(prob,GCE2(Λ=1),tspan;tsargs...);

sols = [sol1,sol2,sol3];

write(prob,eqs[1],sols[1],dn="validate/stochastic/",fn=labels[1])
write(prob,eqs[2],sols[2],dn="validate/stochastic/",fn=labels[2])
write(prob,eqs[3],sols[3],dn="validate/stochastic/",fn=labels[3])
