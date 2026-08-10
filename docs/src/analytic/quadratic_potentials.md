# Quadratic Potentials

```@docs
QuadraticPotentials
```

## Constructing the Field

```@example quadratic
using CairoMakie
using ElectromagneticFields

equ = QuadraticPotentials.init()
```

## Plotting

The three components of the vector potential in the plane ``z = 0``. The first two are linear in
the coordinates, the third is quadratic:

```@example quadratic
plot_equilibrium(equ)
```

## Evaluating the Field

This is one of the few fields that carries an electrostatic potential as well as a magnetic one,
so the generated code includes `φ` and the components of the electric field:

```@example quadratic
QuadraticPotentials.@code()
nothing # hide
```

```@example quadratic
t = 0.0
x = [0.5, 0.3, 0.2]

φ(t, x), [E₁(t, x), E₂(t, x), E₃(t, x)]
```
