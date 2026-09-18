# Perturbations

A perturbation is a field defined in the same way as an equilibrium — by whichever of `A₁`, `A₂`,
`A₃` and `φ` it contributes — and added to one before any code is generated. The sum is traced and
differentiated as a whole, so a perturbed field is a [`FieldFunctions`](@ref) like any other and
costs nothing extra to evaluate.

A perturbation must live on the same chart as the equilibrium it perturbs. The generator asserts
that the two agree on the volume element and the metric, so combining a cartesian perturbation with
a toroidal equilibrium fails immediately rather than producing a field that is quietly wrong.

The package ships one perturbation; writing another is described in
[Adding a Field](@ref).

```@docs
EzCosZPerturbation
```

## Using a Perturbation

The perturbation is the second, optional argument to `FieldFunctions`:

```@example perturbation
using ElectromagneticFields

equ = ThetaPinchEquilibrium()
pert = EzCosZPerturbation(2.0)

plain = FieldFunctions(equ)
perturbed = FieldFunctions(equ, pert)
nothing # hide
```

`EzCosZPerturbation` contributes only a scalar potential, so the magnetic field is untouched:

```@example perturbation
t = 0.0
ξ = [0.5, 0.5, 0.25]

unchanged = B♭(perturbed, t, ξ) ≈ B♭(plain, t, ξ)
@assert unchanged # hide
unchanged
```

while the electrostatic potential and the electric field are not. The θ-pinch alone has neither:

```@example perturbation
φ(plain, t, ξ), φ(perturbed, t, ξ)
```

```@example perturbation
[E♭(plain, t, ξ) E♭(perturbed, t, ξ)]
```

The objects the field was built from remain available, which is how a script can record what it
ran:

```@example perturbation
equilibrium(perturbed), perturbation(perturbed)
```
