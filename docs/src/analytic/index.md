# Analytic Fields

ElectromagneticFields.jl provides a collection of analytically known electromagnetic fields.
Each of them lives in its own submodule, is constructed by an `init` function taking the
parameters of the field, and generates its evaluation routines through a `@code` macro — see
[Usage](../usage.md) for the general workflow, [Coordinates](../coordinates.md) for the charts
listed below, [Fields](../fields.md) for the field quantities each of them provides, and
[Plotting](../plotting.md) for how the figures on the following pages are made.

| Field | Coordinates | |
|---|---|---|
| [Arnold-Beltrami-Childress Field](@ref) | ``(x,y,z)`` | a three-dimensional chaotic field |
| [Axisymmetric Tokamak (Cartesian)](@ref) | ``(x,y,z)`` | simple tokamak with circular flux surfaces |
| [Axisymmetric Tokamak (Cylindrical)](@ref) | ``(R,Z,\phi)`` | the same field in cylindrical coordinates |
| [Axisymmetric Tokamak (Toroidal)](@ref) | ``(r,\theta,\phi)`` | the same field in toroidal coordinates |
| [Dipole](@ref) | ``(x,y,z)`` | magnetic dipole |
| [Penning Traps](@ref) | ``(x,y,z)`` | three Penning trap configurations |
| [Quadratic Potentials](@ref) | ``(x,y,z)`` | electromagnetic field with quadratic potentials |
| [Singular Field](@ref) | ``(x,y,z)`` | magnetic field diverging on the axis |
| [Solov'ev Equilibrium](@ref) | ``(R,Z,\phi)`` | flexible Grad-Shafranov solutions, with and without X-point |
| [Symmetric Solov'ev Equilibrium](@ref) | ``(x,y,z)`` | up/down and left/right symmetric Solov'ev solution |
| [Symmetric Quadratic Field](@ref) | ``(x,y,z)`` | magnetic field growing quadratically with the radius |
| [Theta Pinch](@ref) | ``(x,y,z)`` | homogeneous field along the ``z`` axis |
