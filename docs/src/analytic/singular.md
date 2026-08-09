# Singular Field

```@docs
Singular
```

## Constructing the Field

```@example singular
using CairoMakie
using ElectromagneticFields

equ = Singular.init()
```

## Plotting

Both the vector potential and the magnetic field diverge like ``r^{-3}`` as the ``z`` axis is
approached, so linearly spaced contour levels would show nothing but the singularity. The plot
therefore uses logarithmically spaced levels:

```@example singular
plot_equilibrium(equ)
```

## Evaluating the Field

```@example singular
Singular.@code()
nothing # hide
```

The divergence is steep enough to be worth keeping in mind when this field is used as a test
case — an order of magnitude closer to the axis means three orders of magnitude in ``|B|``:

```@example singular
[B(0.0, r, 0.0, 0.0) for r in (1.0, 0.1, 0.01)]
```
