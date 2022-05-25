using ZonalFlow

tspan   = (0.0,500.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=10000,
            save_everystep=false,
            saveat=50,
            save_noise=false
           );

domain  = Domain(extent=(2π,2π),res=(12,12));
coeffs  = Coefficients(Ω=5.0,θ=0.0,μ=0.01,ν=0.0,ν₄=1.0);
forcing = Stochastic(kf=5,dk=1,ε=0.01);
prob    = BetaPlane(domain,coeffs,forcing);

eqs = [NL(),GQL(0),CE2(),GQL(1),GCE2(1)];
labels = ["nl_sf","ql_sf","ce2_sf","gql_sf","gce2_sf"];

for (eq,label) in zip(eqs,labels)
    sol = integrate(prob,eq,tspan;tsargs...);
    write(prob,eq,sol,dn="data/",fn=label)
end
