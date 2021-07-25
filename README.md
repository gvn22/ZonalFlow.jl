# ZonalFlow
[![Build Status](https://travis-ci.com/gvn22/ZonalFlow.jl.svg?branch=master)](https://travis-ci.com/gvn22/ZonalFlow.jl)
[![Coverage](https://codecov.io/gh/gvn22/ZonalFlow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gvn22/ZonalFlow.jl)

This is a numerical solver for barotropic vorticity equation on the beta-plane, which can solve tendencies for the following dynamical equation systems:

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

### Installation instructions
Simply add the package using the Julia package manager as:

```julia
julia>]
(v1.5) pkg> add ZonalFlow
```

### Example simulation script

```
using ZonalFlow

tspan   = (0.0,1000.0);
tsargs  = (
            dt=0.001,
            adaptive=false,
            progress=true,
            progress_steps=10000,
            save_everystep=false,
            save_start=true,
            dense=false,
            save_noise=false,
            saveat=5
           );

domain  = Domain(extent=(2π,2π),res=(16,16));
coeffs  = Coefficients(Ω=2π,θ=0.0,μ=0.05,ν=0.0,ν₄=0.0);
forcing = PointJet(Ξ=1.0,Δθ=0.1,τ=20.0);
prob    = BetaPlane(domain,coeffs,forcing);

dn      = "pointjet/"

sol     = integrate(prob,NL(),tspan;tsargs...)
write(prob,NL(),sol,dn=dn,fn="nl")


```

### Citing us
We will soon upload citation instructions, as the relevant articles are due for submission soon. Please feel free to contact me if you are using the code!

### License
MIT
