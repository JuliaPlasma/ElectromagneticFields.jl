# Interface

This page is the reference for everything a [`FieldFunctions`](@ref) object provides: the naming
scheme, the complete list of quantities, the shapes they return, the data they carry alongside,
and how to hand a field to code that consumes it.

```@example interface
using ElectromagneticFields
using LinearAlgebra

field = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 2.0))

t = 0.0
ξ = [6.5, 0.5, 0.25]
nothing # hide
```


## The Naming Scheme

A quantity that exists in more than one representation carries the musical isomorphism that
produces that representation. There are three marks, all single Unicode codepoints, all
tab-completable in the REPL and in the editors that support LaTeX completion:

| mark | completion | meaning |
|:--|:--|:--|
| `♭` | `\flat` | index lowered — the covariant components, a one-form |
| `♯` | `\sharp` | index raised — the contravariant components, a vector |
| `♮` | `\natural` | neither — the components in the physical, orthonormal frame |

The bare letter is reserved for the magnitude, which has no index to raise or lower. So `B` is
``|B|``, a scalar, while `B♭`, `B♯` and `B♮` are the three representations of the vector. Two
lowered indices are written with two flats: `B♭♭` is the magnetic two-form ``B_{ij}``.

A leading `D` is a derivative with respect to the chart coordinates, and it is applied to the
quantity the rest of the name denotes. `DA♭` is ``\partial_j A_i``, a matrix; `DDA♭` is
``\partial_k \partial_j A_i``, a rank-three array; `DB` and `DDB` are the gradient and Hessian of
the scalar `B`.

!!! note
    Differentiation is always with respect to the chart's own coordinates ``\xi``, never the
    cartesian ones.


## The Quantities

Every generic takes the field as its first argument, then the time, then the coordinates, given
either as a vector or as three scalars:

```@example interface
B♭(field, t, ξ) == B♭(field, t, ξ...)
```

The return shape is fixed by the tensor rank. Scalars come back as `Float64` (more precisely, as
`float(eltype(ξ))`), and everything else as a `StaticArray`, so extracting a component is an
ordinary index and costs nothing:

```@example interface
B♭(field, t, ξ)[3]
```

### Potentials and fields

| name | shape | quantity |
|:--|:--|:--|
| `A♭`, `A♯` | `SVector{3}` | vector potential, covariant and contravariant |
| `DA♭` | `SMatrix{3,3}` | ``\partial_j A_i`` |
| `DDA♭` | `SArray{3,3,3}` | ``\partial_k \partial_j A_i`` |
| `φ` | scalar | electrostatic potential |
| `B` | scalar | ``\|B\|`` |
| `DB`, `DDB` | `SVector{3}`, `SMatrix{3,3}` | gradient and Hessian of ``\|B\|`` |
| `B♭`, `B♯`, `B♮` | `SVector{3}` | magnetic field, three representations |
| `B♭♭` | `SMatrix{3,3}` | magnetic two-form ``B_{ij}`` |
| `DB♭` | `SMatrix{3,3}` | ``\partial_j B_i`` |
| `b♭`, `b♯`, `b♮` | `SVector{3}` | unit vector along ``B`` |
| `Db♭`, `Db♮` | `SMatrix{3,3}` | ``\partial_j b_i`` and its physical counterpart |
| `DDb♭` | `SArray{3,3,3}` | ``\partial_k \partial_j b_i`` |
| `a♭`, `a♯`, `a♮` | `SVector{3}` | first vector perpendicular to ``b`` |
| `c♭`, `c♯`, `c♮` | `SVector{3}` | second vector perpendicular to ``b`` |
| `E♭`, `E♯` | `SVector{3}` | electric field |
| `DE♭` | `SMatrix{3,3}` | ``\partial_j E_i`` |

### Chart and metric

| name | shape | quantity |
|:--|:--|:--|
| `g♭`, `g♯` | `SMatrix{3,3}` | metric and its inverse |
| `Dg♭`, `Dg♯` | `SArray{3,3,3}` | ``\partial_k g_{ij}`` and ``\partial_k g^{ij}`` |
| `DDg♭`, `DDg♯` | `SArray{3,3,3,3}` | their second derivatives |
| `DF`, `DF̄` | `SMatrix{3,3}` | tangent map ``\partial x^i / \partial \xi^j`` and its inverse |
| `J` | scalar | volume element ``\sqrt{\|g\|} = \|\det DF\|`` |
| `to_cartesian`, `from_cartesian` | `SVector{3}` | the chart map and its inverse |
| `rangemin`, `rangemax` | `SVector{3}` | bounds of the coordinate domain |

`rangemin` and `rangemax` also take a call without the time, since no chart in this package has a
time-dependent domain:

```@example interface
rangemin(field, ξ) == rangemin(field, t, ξ)
```

`functions(field)` returns all of the above at once, as a `NamedTuple` keyed by these names. It is
the same interface `GeometricEquations` uses for the functions of an equation.


## Data Carried Alongside

Four things are values rather than functions, because they do not depend on ``t`` or ``\xi``.

**`parameters(field)`** gives the equilibrium's parameters. The generated code does *not* have
these baked in — it takes them as an argument, and the field supplies them from the equilibrium
struct on every call, which is why two fields of the same type share their compiled code. This is
the place to read a parameter back when a script needs it:

```@example interface
parameters(field)
```

**`coordinates(field)`** gives the equilibrium's own coordinate helpers. Which names are present
varies from field to field — the cartesian fields define `X`, `Y` and `Z`, the axisymmetric ones
add `R`, `r`, `θ` and `ϕ` — which is why they are a `NamedTuple` rather than generics of their
own. Each entry is called like any other field function, minus the field argument:

```@example interface
keys(coordinates(field))
```

```@example interface
coordinates(field).R(t, ξ)
```

**`orientation(field)`** is `+1` or `-1`, the handedness of the chart. See
[Volume Element and Orientation](@ref).

**`periodic(field)`** is an `SVector{3, Bool}`, one entry per coordinate of the chart, saying
whether that coordinate is periodic. It is a property of the chart, so it is answered per chart
family and not per equilibrium, and there is no default — a chart nobody has answered for raises a
`MethodError` when a field is built from it.

This field is on the ``(R, Z, \phi)`` chart, where the toroidal angle is periodic and the two
poloidal coordinates are unbounded:

```@example interface
periodic(field)
```

It is separate from `rangemin`/`rangemax` because a bounded range does not imply periodicity. The
two happen to agree on every chart this package ships — each one bounds exactly its angles, to
``[0, 2\pi]`` — but a chart with a wall at ``r = a``, or a slab bounded in ``z``, would have a
finite range in a coordinate that does not wrap. A chart of your own must therefore state
periodicity itself rather than have it read off the bounds.

The name is `GeometricBase.periodic`, which already means one `Bool` per component. It is
deliberately not `periodicity`: `GeometricEquations` answers that generic with an `(xmin, xmax)`
tuple naming the periodic domain, and a `Bool` vector under the same name would give one generic
two shapes. A `GeometricEquations` problem still takes its own `periodicity` in its own form;
`periodic(field)` is what tells you which components belong in it.

**`equilibrium(field)`** and **`perturbation(field)`** return the objects the field was built
from.


## Handing a Field to a Solver

A `FieldFunctions` is an ordinary immutable value, and every accessor is type-stable and
allocation-free, so it can be stored in the `parameters` of a `GeometricEquations` problem and
used from inside the vector field. This is the intended way to write an equation that depends on
a field:

```julia
using GeometricEquations

function guiding_centre_v(v, t, q, p, params)
    B̄ = B(params.field, t, q)
    b = b♯(params.field, t, q)
    ...
end

prob = ODEProblem(guiding_centre_v, timespan, timestep, q₀;
                  parameters = (field = field, μ = 1.0))
```

Nothing is captured in a closure and nothing is looked up in a global, so the same vector field
works for every equilibrium, and switching fields is a change of `parameters` rather than a
recompilation of the problem.

Two properties make this practical:

* **No allocations.** Every generic returns a `StaticArray` or a scalar, whatever container the
  coordinates arrive in.
* **Type stability.** Every quantity of one field returns the same element type, `float` of the
  coordinates' own type, whether or not its expression happens to mention the coordinates.

For every equilibrium the test suite measures the first, and for the second compares the element
type of each value returned. Neither bullet is asserted through inference: `@inferred` appears
nowhere in the suite.

Measuring the first is worth a word of warning: `@allocated` charges for the lookup of a
non-`const` global, so it has to be done behind a function barrier or it reports the boxing of
`field` rather than the work of the call.

```@example interface
function measure(f, field, t, ξ)
    f(field, t, ξ)                      # warm up
    @allocated f(field, t, ξ)
end

measure(DDg♯, field, t, ξ), eltype(g♭(field, t, ξ))
```


## Identities

The relations between the representations hold pointwise, which is worth knowing both as
documentation and as the check to run when adding a chart:

```@example interface
G, Ḡ = g♭(field, t, ξ), g♯(field, t, ξ)
F, F̄ = DF(field, t, ξ), DF̄(field, t, ξ)

identities = (
    Ḡ ≈ inv(G),                                     # the metric and its inverse
    F̄ ≈ inv(F),                                     # the tangent map and its inverse
    F' * F ≈ G,                                     # the metric is the pullback
    J(field, t, ξ) ≈ sqrt(det(F' * F)),             # the volume element
    det(F) ≈ orientation(field) * J(field, t, ξ),   # the signed determinant
    B♯(field, t, ξ) ≈ Ḡ * B♭(field, t, ξ),          # raising an index
    B♮(field, t, ξ) ≈ F̄' * B♭(field, t, ξ),         # the physical components
    B(field, t, ξ) ≈ sqrt(B♯(field, t, ξ)' * B♭(field, t, ξ)),
    b♭(field, t, ξ) ≈ B♭(field, t, ξ) / B(field, t, ξ),
    B♭♭(field, t, ξ) ≈ -B♭♭(field, t, ξ)'           # the two-form is antisymmetric
)

@assert all(identities) # hide
all(identities)
```

The perpendicular frame is orthonormal, which is visible directly in the physical components:

```@example interface
frame = [a♮(field, t, ξ) b♮(field, t, ξ) c♮(field, t, ξ)]

@assert frame' * frame ≈ I # hide
round.(frame' * frame; digits = 12)
```
