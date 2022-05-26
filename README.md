# ZonalFlow
[![build](https://github.com/gvn22/ZonalFlow.jl/actions/workflows/ci.yml/badge.svg?event=push)](https://github.com/gvn22/ZonalFlow.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/gvn22/ZonalFlow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gvn22/ZonalFlow.jl)

ZonalFlow is a spectral solver for barotropic vorticity equation on the beta-plane, which can solve tendencies for the following dynamical equation systems:

- NL: original (fully non-linear) governing equations with no mode-reduction
- QL: quasilinear equations
- GQL: generalised quasilinear equations

and the following (equivalent) statistical equation systems:

- CE2: cumulant expansions at second order
- GCE2: generalized cumulant expansions at second order

The package interfaces with the DifferentialEquations package in order to utilise its ecosystem of time integration schemes.  A variety of different initial condition types can  be chosen. Simulation data is output as npz files, which can be accessed using Python.

## Contents
* [Installation instructions](#installation-instructions)
* [Example simulation script](#example-simulation-script)
* [Citing us](#citing-us)
* [License](#license)

### Installation
Add ZonalFlow using the Julia package manager as:

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
