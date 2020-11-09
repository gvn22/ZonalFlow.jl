# ZonalFlow
[![Build Status](https://travis-ci.com/gvn22/ZonalFlow.jl.svg?branch=master)](https://travis-ci.com/gvn22/ZonalFlow.jl)
[![Coverage](https://codecov.io/gh/gvn22/ZonalFlow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gvn22/ZonalFlow.jl)

## Generalized Quasilinear approximation and Generalized Cumulant Expansion

- `nl()`: Numerically solving the fully non-linear (NL) dynamics of dissipative and driven rotational flows can be a computational expensive task if one is interested in the statistical behaviour of turbulent zonal jets. In a spectral code, the computational bottleneck is computing the sum of all interactions arising from the non-linear terms in the dynamical field equations.

- `gql()`: The Generalized Quasilinear (GQL) approximation [[1]](#1) simplifies these non-linear terms by inducing rules for interaction between low (L) and high (H) wavenumber zonal projections of the dynamical field obtained using a spectral filter with cutoff Λ. Specifically, the `HH→H` (eddy-eddy non-linearity interactions as well as the `HL→L` and `LL→H` interactions are eliminated. In this manner, the field equations can be extrapolated between Quasilinear (QL) dynamics for Λ = 0, and NL dynamics for Λ = M (the maximum zonal wavenumber). Therefore, a GQL system obtained for Λ < M, and that suffices to simulate zonal jet statistics for given friction and driving parameters, is by nature a reduced model of the underlying dynamics. However, obtaining statistics of the flow may still require performing simulations with large spin-up times to arrive at a statistical significant sample.

- `gce2()`: The Generalized Cumulant Expansion (GCE2) circumvents this last problem, by posing equations in terms of the statistics -- the first and second cumulant -- derived from the GQL equations for a given cutoff Λ. This allows the required cumulants (or equaivalently the moments) to be obtained directly, precluding the need for large spin-up dynamical calculations. In practice, the cumulant sizes for low Λ are large and can be computationally more expensive per timestep; however, this can be partly offset by the fact that fewer timesteps need to be solved. Use of dimensional reduction techniques can help speed this up further.

`ZonalFlow` allows you to solve all three sets of equations: NL, GQL and GCE2 for dynamics on the β-plane. It uses the [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl) package for numerically time-integrating the spectral ODE problem, giving access to a range of integration algorithms. A matrix-free representation is used to minimize computational cost. Currently, only fixed timestep algorithms are recommended and computations are serial for now.

Functions `nl(params...)`, `gql(params...)` and `gce2(params...)` solver functions are exported by the package together with a number of solution analysis functions.

###### References
<a id="1">[1]</a> Marston, J. B. and Chini, G. P. and Tobias, S. M. (2016) _Physical Review Letters_ __116__ 214501

## Use
The current release solves for equations on the β-plane with forcing by a deterministic point-jet and relaxation; the choice of hyperviscosity is available. Set the parameters for a given solution as follows:

Domain:
- lx: axial domain length
- ly: transverse domain length
- nx: axial resolution
- ny: transverse resolution
- Λ (for GQL and GCE2)

Problem:
- Ω: Rotational rate
- θ: Latitude
- β: Coriolis parameter
- Ξ: jet strength
- τ: relaxation time

Timestepping and solution save:
- jw: Jet width; defaults to 0.2
- ic: Initial Condition; defaults to point jet with random seed.
- dt: Timestep size; defaults to 0.001.
- t_end: Final solution time
- savefreq: Data saving frequency

Each solver function returns the solution variable to which the following analysis algorithms may be applied.

- `zonalenergy(params...)`
- `fourierenergy(params...)`
- `meanvorticity(params...)`
- `zonalvelocity(params...)`

## License
MIT
