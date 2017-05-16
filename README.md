
# MagneticEquilibria.jl

*Common Interface for Magnetic Equilibria*

MagneticEquilibria.jl provides a common interface for evaluating analytical and numerical equilibria.
For analytical equilibria, it generates Julia code providing high-level evaluation routines. Numerical equilibria
are interpolated using ApproXD and can be evaluated with the very same interface as the analytical equilibria.


## Features

Analytical equilibria:
- simple axisymmetric tokamak equilibrium
- flexible Solov'ev equilibria with and without X-point
- 3D perturbations (e.g., magnetic islands; planned)

Numerical equilibria (planned):
- projected analytic equilibria
- EFIT
- VMEC


## License

The MagneticEquilibria.jl package is licensed under the [MIT "Expat" License](LICENSE.md).


## References

- Antoine J. Cerfon, Jeffrey P. Freidberg. "One size fits all" analytic solutions to the Grad–Shafranov equation. Physics of Plasmas 17 (3), 032502.
