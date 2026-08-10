# Theta Pinch

```@docs
ThetaPinch
```

## Constructing the Field

```@example thetapinch
using CairoMakie
using ElectromagneticFields

equ = ThetaPinch.init()
```

## Plotting

The magnetic field is homogeneous and points along ``z``, so the two components of the vector
potential are linear in ``y`` and ``x`` respectively:

```@example thetapinch
plot_equilibrium(equ)
```

This is the simplest field in the collection and therefore a good starting point when checking
that a new integrator or diagnostic behaves as expected.
