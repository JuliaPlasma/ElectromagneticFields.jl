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

The flux surfaces of this equilibrium are circular. In the ``(x,z)`` plane at ``y = 0`` they show
up as the contours of the ``y`` component of the vector potential:

```@example atc
plot_equilibrium(equ)
```

The same field in [cylindrical](@ref "Axisymmetric Tokamak (Cylindrical)") and
[toroidal](@ref "Axisymmetric Tokamak (Toroidal)") coordinates is available as well, which makes
this a convenient test case for coordinate transformations.
