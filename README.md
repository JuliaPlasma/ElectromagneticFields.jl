
# ElectromagneticFields.jl

*Common Interface for Electromagnetic Fields*

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaplasma.github.io/ElectromagneticFields.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliaplasma.github.io/ElectromagneticFields.jl/latest/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/E/ElectromagneticFields.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/E/ElectromagneticFields.html)
![CI](https://github.com/JuliaPlasma/ElectromagneticFields.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/JuliaPlasma/ElectromagneticFields.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaPlasma/ElectromagneticFields.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3662494.svg)](https://doi.org/10.5281/zenodo.3662494)

ElectromagneticFields.jl provides a common interface for evaluating analytical and numerical magnetic equilibria, general electromagnetic fields and in the future also simple Maxwell solvers.
For analytical equilibria, it generates Julia code providing high-level evaluation routines. Numerical equilibria
are interpolated using ApproXD and can be evaluated with the very same interface as the analytical equilibria.


## Features

Analytical equilibria:
- simple axisymmetric tokamak equilibrium in cartesian, cylindrical and toroidal coordinates
- flexible Solov'ev equilibria with and without X-point (including ITER, NSTX and FRC configurations)
- symmetric quadratic and singular magnetic fields, symmetric Solov'ev equilibria
- Penning trap
- Arnold-Beltrami-Childress (ABC) field
- 3D perturbations (e.g., magnetic islands, electric fields)

Numerical equilibria (planned):
- projected analytic equilibria
- EFIT
- VMEC

Numerical solvers (planned):
- B-Spline, FEM and pseudo-spectral Poisson, Ampère and Faraday solvers


## References

- Antoine J. Cerfon, Jeffrey P. Freidberg. "One size fits all" analytic solutions to the Grad–Shafranov equation. [Physics of Plasmas 17, 032502, 2010](https://doi.org/10.1063/1.3328818).
- Patrick J. McCarthy. Analytical solutions to the Grad–Shafranov equation for tokamak equilibrium with dissimilar source functions. [Physics of Plasmas 6, 3554, 1999](https://doi.org/10.1063/1.873630).
- Yanyan Shi, Yajuan Sun, Yulei Wang, Jian Liu, Study of adaptive symplectic methods for simulating charged particle dynamics, [Journal of Computational Dynamics 6, 429-448, 2019](http://dx.doi.org/10.3934/jcd.2019022).


## License

The ElectromagneticFields.jl package is licensed under the [MIT "Expat" License](LICENSE.md).
