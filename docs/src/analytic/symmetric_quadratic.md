# Symmetric Quadratic Field

```@docs
SymmetricQuadratic
```

## Constructing the Field

```@example symmetric
using CairoMakie
using ElectromagneticFields

equ = SymmetricQuadratic.init(1.0)
```

## Plotting

The two components of the vector potential and the ``z`` component of the magnetic field, which
grows quadratically with the distance from the axis:

```@example symmetric
plot_equilibrium(equ)
```

## Evaluating the Field

```@example symmetric
SymmetricQuadratic.@code(1.0)
nothing # hide
```

The field is axisymmetric, so ``|B|`` depends on the radius alone:

```@example symmetric
B(0.0, 0.3, 0.4, 0.0) ≈ B(0.0, 0.5, 0.0, 0.0)
```
