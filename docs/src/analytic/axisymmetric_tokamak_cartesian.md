# Axisymmetric Tokamak (Cartesian)

```@docs
AxisymmetricTokamakCartesian
```

## Constructing the Field

```@example atc
using CairoMakie
using ElectromagneticFields

equ = AxisymmetricTokamakCartesian.init()
```

## Plotting

The flux surfaces of this equilibrium are circular. In the ``(x,z)`` plane at ``y = 0`` the ``y``
component of the vector potential is the physical toroidal component ``A_\phi``, so ``R \, A_y`` is
the poloidal flux function and its contours are the flux surfaces:

```@example atc
plot_equilibrium(equ)
```

The same field in [cylindrical](@ref "Axisymmetric Tokamak (Cylindrical)") and
[toroidal](@ref "Axisymmetric Tokamak (Toroidal)") coordinates is available as well, which makes
this a convenient test case for coordinate transformations.
