# Axisymmetric Tokamak (Toroidal)

```@docs
AxisymmetricTokamakToroidalEquilibrium
```

## Constructing the Field

```@example att
using CairoMakie
using ElectromagneticFields

equ = AxisymmetricTokamakToroidalEquilibrium()
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
field = FieldFunctions(AxisymmetricTokamakToroidalEquilibrium())
nothing # hide
```

The chart maps ``(r, \theta, \phi)`` to the cartesian coordinates, and `to_cartesian` and
`from_cartesian` move between the two:

```@example att
t = 0.0
ξ = [0.5, π/4, 0.0]

x = to_cartesian(field, t, ξ)
```

```@example att
roundtrip = from_cartesian(field, t, x) ≈ ξ
@assert roundtrip # hide
roundtrip
```

Because the toroidal chart is left-handed, the determinant of its Jacobian is minus the Jacobian
determinant `J`:

```@example att
using LinearAlgebra

signed = det(DF(field, t, ξ)) ≈ orientation(field) * J(field, t, ξ)
@assert signed # hide
orientation(field), signed
```


## The Regularized Gauge

```@docs
AxisymmetricTokamakToroidalRegularizationEquilibrium
```

The vector potential is fixed only up to a gauge, and the two gauges are genuinely different
equilibria here, because the generated code follows whatever was written down. The poloidal
vector potential of the chart above is singular on the magnetic axis; this variant carries the
same magnetic field in a gauge that is regular there.

```@example att
regular = FieldFunctions(AxisymmetricTokamakToroidalRegularizationEquilibrium())

same_field = B♭(regular, t, ξ) ≈ B♭(field, t, ξ)
@assert same_field # hide
same_field
```

The vector potentials differ, as they must:

```@example att
[A♭(field, t, ξ) A♭(regular, t, ξ)]
```
