# Dipole

```@docs
Dipole
```

## Constructing the Field

The only parameter is the field strength ``B_0``, which defaults to `1000.0`:

```@example dipole
using CairoMakie
using ElectromagneticFields

equ = Dipole.init()
```

## Plotting

The vector potential circles the ``z`` axis and falls off with the third power of the distance
from the origin. The plot shows its two non-vanishing components in the plane ``z = 1``:

```@example dipole
plot_equilibrium(equ)
```
