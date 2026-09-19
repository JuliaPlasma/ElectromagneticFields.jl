# ElectromagneticFields.jl

*Common Interface for Electromagnetic Fields*

[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/E/ElectromagneticFields.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/E/ElectromagneticFields.html)
[![CI](https://github.com/JuliaPlasma/ElectromagneticFields.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaPlasma/ElectromagneticFields.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/JuliaPlasma/ElectromagneticFields.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/JuliaPlasma/ElectromagneticFields.jl/actions/workflows/Documentation.yml)
[![Coverage](https://codecov.io/gh/JuliaPlasma/ElectromagneticFields.jl/graph/badge.svg?token=shiEHXD1rj)](https://codecov.io/gh/JuliaPlasma/ElectromagneticFields.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3662494.svg)](https://doi.org/10.5281/zenodo.3662494)

ElectromagneticFields.jl provides a collection of analytically known electromagnetic fields —
tokamak equilibria, Solov'ev solutions of the Grad-Shafranov equation, Penning traps, the
Arnold-Beltrami-Childress field and several others — together with the machinery that turns any
of them into fast, allocation-free evaluation code.

A field is described by very little: the covariant components of its vector potential, an
electrostatic potential where there is one, and the chart it is written in. Everything a
simulation actually needs — the magnetic field, its magnitude, the unit vector along it, an
orthonormal frame adapted to it, the electric field, the metric, and the derivatives of all of
these — is derived from that symbolically and compiled.

```@example index
using ElectromagneticFields

field = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 2.0))

B♭(field, 0.0, [6.5, 0.5, 0.25])
```

The result is a value. It can be built inside a function, stored, and passed to a vector field
that depends on it — typically in the `parameters` of a `GeometricEquations` problem, from which
the equation function evaluates the components it needs.

Where to go from here:

* [Usage](usage.md) — constructing an equilibrium and generating its code.
* [Interface](interface.md) — the complete list of quantities, their shapes, and how to hand a
  field to a solver.
* [Code Generation](generation.md) — how the code is derived, what it costs, and how to add a
  field of your own.
* [Coordinates](coordinates.md) and [Fields](fields.md) — the differential geometry behind the
  naming and the derivation.
* [Analytic Fields](analytic/index.md) — the fields the package ships.


## References

- Antoine J. Cerfon, Jeffrey P. Freidberg. "One size fits all" analytic solutions to the Grad–Shafranov equation. [Physics of Plasmas 17, 032502, 2010](https://doi.org/10.1063/1.3328818).
- Patrick J. McCarthy. Analytical solutions to the Grad–Shafranov equation for tokamak equilibrium with dissimilar source functions. [Physics of Plasmas 6, 3554, 1999](https://doi.org/10.1063/1.873630).
- Yanyan Shi, Yajuan Sun, Yulei Wang, Jian Liu, Study of adaptive symplectic methods for simulating charged particle dynamics, [Journal of Computational Dynamics 6, 429-448, 2019](http://dx.doi.org/10.3934/jcd.2019022).


## License

> Copyright (c) Michael Kraus <michael.kraus@ipp.mpg.de>
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.
