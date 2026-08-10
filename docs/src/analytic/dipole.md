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

The vector potential circles the ``z`` axis with a magnitude ``B_0 \, \rho / r^3``, where ``\rho``
is the distance from the axis and ``r`` the distance from the origin. The plot shows its two
non-vanishing components in the plane ``z = 1``:

```@example dipole
plot_equilibrium(equ)
```
