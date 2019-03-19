
# ElectromagneticFields.jl

*Common Interface for Electromagnetic Fields*

ElectromagneticFields.jl provides a common interface for evaluating analytical and numerical magnetic equilibria, general electromagnetic fields and in the future also simple Maxwell solvers.
For analytical equilibria, it generates Julia code providing high-level evaluation routines. Numerical equilibria
are interpolated using ApproXD and can be evaluated with the very same interface as the analytical equilibria.


## Features

Analytical equilibria:
- simple axisymmetric tokamak equilibrium in cartesian, cylindrical and toroidal coordinates
- flexible Solov'ev equilibria with and without X-point
- simple quadratic Solov'ev equilibrium as well as symmetric quadratic field
- ABC field
- 3D perturbations (e.g., magnetic islands)

Numerical equilibria (planned):
- projected analytic equilibria
- EFIT
- VMEC

Numerical solvers (planned):
- B-Spline, FEM and pseudo-spectral Poisson, Ampère and Faraday solvers


## License

The ElectromagneticFields.jl package is licensed under the [MIT "Expat" License](LICENSE.md).


## References

- Antoine J. Cerfon, Jeffrey P. Freidberg. "One size fits all" analytic solutions to the Grad–Shafranov equation. Physics of Plasmas 17 (3), 032502.
