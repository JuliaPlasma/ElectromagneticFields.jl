# Usage

Working with ElectromagneticFields.jl has two steps. First an equilibrium is constructed, which
is a small immutable object holding nothing but the parameters of the field. Then the code for
evaluating that field is generated from it into a [`FieldFunctions`](@ref) object, and that is
what the actual computation uses.


## Constructing an Equilibrium

Every field is a type with sensible defaults, so the shortest way to get an equilibrium is

```@example usage
using ElectromagneticFields

equ = AxisymmetricTokamakCylindricalEquilibrium()
```

All parameters can be passed explicitly, here the major radius ``R_0``, the magnetic field
strength ``B_0`` at the magnetic axis, and the safety factor ``q_0``:

```@example usage
equ = AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 2.0)
```

Some fields ship named configurations as well, e.g. the Solov'ev equilibrium comes with parameter
sets for ITER, NSTX and a field reversed configuration:

```@example usage
SolovevEquilibriumITER()
```


## Generating the Evaluation Code

An equilibrium object on its own cannot be evaluated. The evaluation routines are traced
symbolically from its chart, metric and vector potential, differentiated where necessary, and
compiled. [`FieldFunctions`](@ref) does all of that and returns the result as a value:

```@example usage
field = FieldFunctions(equ)
nothing # hide
```

This is an ordinary object. It can be built inside a function, stored in a `Dict`, or passed to a
vector field that depends on the field — typically through a `GeometricEquations` problem's
`parameters`, from which the equation function evaluates the components it needs.

Every quantity is reached through a generic that takes the field first. They all take the time as
their second argument, followed by the coordinates of the evaluation point, which for this
equilibrium are ``(R, Z, \phi)``:

```@example usage
t = 0.0
x = [6.5, 0.5, 0.0]

B(field, t, x)
```

Both a vector and a splatted call are defined, so `B(field, t, x)` and `B(field, t, x...)` are the
same thing:

```@example usage
same = B(field, t, x...) == B(field, t, x)
@assert same # hide
same
```

A perturbation is combined with the equilibrium before any code is generated, so a perturbed field
is a `FieldFunctions` like any other:

```@example usage
perturbed = FieldFunctions(ThetaPinchEquilibrium(), EzCosZPerturbation())
φ(perturbed, t, [0.5, 0.5, 0.25])
```


## What Gets Generated

Rather more than just the magnetic field, and each quantity comes as a whole tensor rather than
one function per component. A vector-valued quantity returns an `SVector{3}`, a Jacobian or the
metric an `SMatrix{3,3}`, and the higher derivatives an `SArray` of the matching rank, so a single
component is an ordinary index:

```@example usage
B♭(field, t, x)[1]
```

The naming follows the musical isomorphisms of differential geometry: `♭` lowers an index and
gives the covariant components, `♯` raises one and gives the contravariant components, and `♮`
gives the components in the physical (orthonormal) frame. They are typed `\flat`, `\sharp` and
`\natural`. The bare letter is the magnitude, so `B` is ``|B|``. See [Coordinates](coordinates.md)
for what the three representations are and [Fields](fields.md) for how the quantities below are
derived from the vector potential.

| | covariant | contravariant | physical | derivatives |
|:--|:--|:--|:--|:--|
| vector potential | `A♭` | `A♯` | | `DA♭`, `DDA♭` |
| magnetic field | `B♭` | `B♯` | `B♮` | `DB♭`; two-form `B♭♭` |
| magnitude of the magnetic field | `B` | | | `DB`, `DDB` |
| unit magnetic field | `b♭` | `b♯` | `b♮` | `Db♭`, `Db♮`, `DDb♭` |
| perpendicular frame | `a♭`, `c♭` | `a♯`, `c♯` | `a♮`, `c♮` | |
| electric field | `E♭` | `E♯` | | `DE♭` |
| metric | `g♭` | `g♯` | | `Dg♭`, `Dg♯`, `DDg♭`, `DDg♯` |

| | |
|---|---|
| `φ` | electrostatic potential |
| `J` | volume element of the chart, ``\sqrt{\|g\|} = \|\det DF\|`` |
| `DF`, `DF̄` | Jacobian matrix of the chart and its inverse |
| `from_cartesian`, `to_cartesian` | coordinate transformations |
| `rangemin`, `rangemax` | bounds of the coordinate domain |

Four things are data rather than functions. The parameters of the equilibrium, which the
generated code has baked in as literals:

```@example usage
parameters(field)
```

the equilibrium's own coordinate helpers, whose names vary from field to field:

```@example usage
keys(coordinates(field))
```

which coordinates of the chart are periodic, `periodic(field)`, and the handedness of the
chart. `J` is the volume element and `DF` the Jacobian matrix; the two are related through
`orientation`, and this equilibrium uses a left-handed ``(R, Z, \phi)`` chart, so the
determinant of `DF` is `-J`:

```@example usage
using LinearAlgebra

signed = det(DF(field, t, x)) ≈ orientation(field) * J(field, t, x)
@assert signed # hide
orientation(field), signed
```

`functions(field)` returns all of the generated functions at once, as a `NamedTuple` keyed by the
names in the tables above.


## Evaluating on a Grid

Since the generated functions are ordinary compiled Julia functions, sampling a field is just a
comprehension:

```@example usage
R₀ = parameters(field).R₀

Rgrid = LinRange(0.5 * R₀, 1.5 * R₀, 100)
Zgrid = LinRange(-0.5 * R₀, 0.5 * R₀, 120)

Bfield = [B(field, t, Rgrid[i], Zgrid[j], 0.0)
          for i in eachindex(Rgrid), j in eachindex(Zgrid)]

size(Bfield)
```

Note the index order: `Bfield[i,j]` holds the value at `(Rgrid[i], Zgrid[j])`, which is also the
convention Makie's `contour` expects, so such an array can be passed straight to a plotting call.
See [Plotting](plotting.md) for how to turn this into a figure.
