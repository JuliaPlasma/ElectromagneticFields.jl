# Penning Traps

A Penning trap confines charged particles with a homogeneous axial magnetic field and a
quadrupole electric field. Three variants are provided, differing in the shape of the magnetic
field: uniform, a magnetic bottle, and an asymmetric configuration.

All three carry an electrostatic potential in addition to the magnetic field, so the generated
code includes `φ` and the components of the electric field.


## Uniform Magnetic Field

```@docs
PenningTrapUniform
```

```@example penning
using ElectromagneticFields

equ = PenningTrapUniform.init()
```

```@example penning
PenningTrapUniform.@code()
nothing # hide
```

The magnetic field is homogeneous along ``z``, while the electric field is linear in all three
coordinates and pulls the particle back towards the mid-plane:

```@example penning
t = 0.0
x = [0.1, 0.2, 0.3]

B(t, x), [E₁(t, x), E₂(t, x), E₃(t, x)]
```


## Magnetic Bottle

```@docs
PenningTrapBottle
```

```@example penning
PenningTrapBottle.init()
```


## Asymmetric Magnetic Field

```@docs
PenningTrapAsymmetric
```

```@example penning
PenningTrapAsymmetric.init()
```


These fields have no dedicated plotting routine. They are simple enough to sample and plot
directly — see [Plotting by Hand](@ref) for how to do that.
