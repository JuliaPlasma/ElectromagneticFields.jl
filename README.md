
# ElectromagneticFields.jl

*Common Interface for Electromagnetic Fields*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaplasma.github.io/ElectromagneticFields.jl/stable/)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliaplasma.github.io/ElectromagneticFields.jl/latest/)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/E/ElectromagneticFields.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/E/ElectromagneticFields.html)
[![CI](https://github.com/JuliaPlasma/ElectromagneticFields.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaPlasma/ElectromagneticFields.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/JuliaPlasma/ElectromagneticFields.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/JuliaPlasma/ElectromagneticFields.jl/actions/workflows/Documentation.yml)
[![Coverage](https://codecov.io/gh/JuliaPlasma/ElectromagneticFields.jl/graph/badge.svg?token=shiEHXD1rj)](https://codecov.io/gh/JuliaPlasma/ElectromagneticFields.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3662494.svg)](https://doi.org/10.5281/zenodo.3662494)

ElectromagneticFields.jl provides a collection of analytical magnetic equilibria and electromagnetic fields.
It generates Julia code providing high-level evaluation routines.


## Features

Analytical equilibria:
- simple axisymmetric tokamak equilibrium in cartesian, cylindrical and toroidal coordinates
- flexible Solov'ev equilibria with and without X-point (including ITER, NSTX and FRC configurations)
- symmetric quadratic and singular magnetic fields, symmetric Solov'ev equilibria
- Penning trap
- Arnold-Beltrami-Childress (ABC) field
- 3D perturbations (e.g., magnetic islands, electric fields)


## References

- Antoine J. Cerfon, Jeffrey P. Freidberg. "One size fits all" analytic solutions to the Grad–Shafranov equation. [Physics of Plasmas 17, 032502, 2010](https://doi.org/10.1063/1.3328818).
- Patrick J. McCarthy. Analytical solutions to the Grad–Shafranov equation for tokamak equilibrium with dissimilar source functions. [Physics of Plasmas 6, 3554, 1999](https://doi.org/10.1063/1.873630).
- Yanyan Shi, Yajuan Sun, Yulei Wang, Jian Liu, Study of adaptive symplectic methods for simulating charged particle dynamics, [Journal of Computational Dynamics 6, 429-448, 2019](http://dx.doi.org/10.3934/jcd.2019022).


## Development

### Git hooks

Two hooks live in `.githooks`. They are **not active in a fresh clone** — `core.hooksPath` is local
configuration and does not travel with a push — so enable them once per clone:

```sh
git config core.hooksPath .githooks
```

**`pre-commit`** acts on **staged `.jl` files only**, and exits immediately when a commit stages
none, so a documentation- or workflow-only commit is not slowed down by it:

- **JuliaFormatter `--check`**, honouring this repository's own `.JuliaFormatter.toml` — **blocks**
  the commit. Formatting is mechanical and always fixable.
- **`fatou lint`**, when `fatou` is installed — **advisory only**, and deliberately so: its
  `unused-import` rule does not follow `include`, so it flags the load-bearing imports of every
  module file.
- **`using <Package>`**, which catches a syntax error or a broken `include` — **blocks**.

**`pre-push`** runs the full test suite with `--check-bounds=auto`, but **only when pushing to
`main` or `master`**; a topic branch is left to CI. It prints nothing for **10–30 minutes**, which
looks exactly like a network hang and is not one. If you do interrupt it, check for an orphaned
Julia process that the killed hook left behind.

Either hook can be bypassed for a single command with `--no-verify`, for a change you know it does
not apply to:

```sh
git commit --no-verify
git push --no-verify
```

The hooks are generated from one shared copy and are byte-identical across the related
repositories, so edit them there rather than here — a local edit is silently undone by the next
install.
