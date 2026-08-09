# Axisymmetric Tokamak (Toroidal)

```@docs
AxisymmetricTokamakToroidal
```

## Constructing the Field

```@example att
using CairoMakie
using ElectromagneticFields

equ = AxisymmetricTokamakToroidal.init()
```

## Plotting

The coordinates ``(r, \theta, \phi)`` are the natural ones for this field — the flux surfaces are
the surfaces of constant ``r``. For the plot the poloidal plane is sampled in ``(R,Z)`` and mapped
to ``(r,\theta)`` first, so the result is directly comparable with the
[cartesian](@ref "Axisymmetric Tokamak (Cartesian)") and
[cylindrical](@ref "Axisymmetric Tokamak (Cylindrical)") versions of the same equilibrium:

```@example att
plot_equilibrium(equ)
```

## Evaluating the Field

```@example att
AxisymmetricTokamakToroidal.@code()
nothing # hide
```

The chart maps ``(r, \theta, \phi)`` to the cartesian coordinates, and `to_cartesian` and
`from_cartesian` move between the two:

```@example att
t = 0.0
ξ = [0.5, π/4, 0.0]

x = to_cartesian(t, ξ)
```

```@example att
from_cartesian(t, x) ≈ ξ
```

Because the toroidal chart is left-handed, the determinant of its Jacobian is minus the Jacobian
determinant `J`:

```@example att
using LinearAlgebra

orientation(), det(DF(t, ξ)) ≈ orientation() * J(t, ξ)
```
