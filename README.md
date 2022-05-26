# ZonalFlow
[![build](https://github.com/gvn22/ZonalFlow.jl/actions/workflows/ci.yml/badge.svg?event=push)](https://github.com/gvn22/ZonalFlow.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/gvn22/ZonalFlow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gvn22/ZonalFlow.jl)

`ZonalFlow` is a spectral solver for Direct Numerical Simulation (DNS) and Direct Statistical Simulation (DSS) of the barotropic vorticity equation on the beta-plane.

DNS can be performed for the following equations:

- NL: original (fully non-linear) master equations
- QL: quasilinear equations
- GQL: generalised quasilinear equations

DSS can be performed making use of the following equations:

- CE2: cumulant expansions at second order
- GCE2: generalised cumulant expansions at second order

`ZonalFlow` uses time integrators from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) and [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) packages. Choice of unity-rank or full-rank initialisation is provided, and the following types of problems (corresponding to different driving mechanisms) can be solved:

* Deterministic pointjet [1](https://journals.ametsoc.org/view/journals/atsc/65/6/2007jas2510.1.xml)
* Two-scale Kolmogorov flow [2](https://aip.scitation.org/doi/10.1063/1.5004683)
* Stochastic narrow-band forcing [3](https://journals.ametsoc.org/view/journals/atsc/73/5/jas-d-15-0288.1.xml)

Simulation data is saved as `.jld2` files and the post-processed output is saved in `.npz` files (sample `Python` post-processing files are located in the examples directory).

## Contents
* [Installation](#installation-instructions)
* [Examples](#example-simulation-script)
* [License](#license)

### Installation
Add `ZonalFlow` using the `Julia` package manager as:

```julia
julia> using Pkg
julia> Pkg.add("ZonalFlow")
```

### Examples

Example scripts are located in the examples directory. A fully-nonlinear solution of stochastically-driven jets can be obtained by the following code:

```
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
```

### License
MIT
