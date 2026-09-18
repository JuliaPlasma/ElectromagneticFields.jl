# Penning Traps

A Penning trap confines charged particles with a homogeneous axial magnetic field and a
quadrupole electric field. Three variants are provided, differing in the shape of the magnetic
field: uniform, a magnetic bottle, and an asymmetric configuration.

All three carry an electrostatic potential in addition to the magnetic field, so the generated
code includes `φ` and the components of the electric field.


## Uniform Magnetic Field

```@docs
PenningTrapUniformEquilibrium
```

```@example penning
using ElectromagneticFields

equ = PenningTrapUniformEquilibrium()
```

```@example penning
field = FieldFunctions(equ)
nothing # hide
```

The magnetic field is homogeneous along ``z``, while the electric field is linear in all three
coordinates and pulls the particle back towards the mid-plane:

```@example penning
t = 0.0
x = [0.1, 0.2, 0.3]

B(field, t, x), E♭(field, t, x)
```


## Magnetic Bottle

```@docs
PenningTrapBottleEquilibrium
```

```@example penning
PenningTrapBottleEquilibrium()
```


## Asymmetric Magnetic Field

```@docs
PenningTrapAsymmetricEquilibrium
```

```@example penning
PenningTrapAsymmetricEquilibrium()
```


These fields have no dedicated plotting routine. They are simple enough to sample and plot
directly — see [Plotting by Hand](@ref) for how to do that.
