# Code Generation

An equilibrium object carries nothing but the parameters of a field. It cannot be evaluated, and
it is not meant to be: the quantities a simulation needs are the magnetic field, the frame along
it and a tower of derivatives, none of which the equilibrium states. They are *derived* from it,
symbolically, and compiled. This page describes how, what that costs, and what to do when the
result is not what was expected.


## What the Equilibrium Supplies

An equilibrium is required to define very little. Everything on this list is a method taking the
coordinates and the equilibrium, written generically enough to accept symbolic arguments:

| method | meaning | default |
|:--|:--|:--|
| `A₁`, `A₂`, `A₃` | covariant components of the vector potential | none, must be defined |
| `φ` | electrostatic potential | `0` |
| `x¹`, `x²`, `x³` | the chart map to cartesian coordinates | none for a new chart |
| `ξ¹`, `ξ²`, `ξ³` | its inverse | none for a new chart |
| `g₁₁` … `g₃₃` | metric coefficients, in the chart's own coordinates | the identity |
| `J` | volume element ``\sqrt{\|g\|}`` | none for a new chart |
| `orientation` | `+1` or `-1` | `+1` |
| `minx¹` … `maxx³` | bounds of the coordinate domain | ``\pm\infty`` |
| `get_functions` | coordinate helpers to expose as `coordinates(field)` | none |
| `get_parameters` | which fields are parameters | all but `name` |

A field built on `CartesianEquilibrium` inherits the identity chart, `J = 1`, the euclidean
metric and `orientation = +1`, so in that case only the vector potential is left. That is why
`ThetaPinchEquilibrium` is a few dozen lines.


## The Symbolic Trace

`FieldFunctions(equ)` begins by declaring symbolic variables for the time and the three
coordinates and calling the methods above with them. Because those methods are written
generically over `AbstractVector`, they return symbolic expressions rather than numbers, and what
comes back is the vector potential and the chart as algebraic expressions in ``\xi``.

From there the derivation is the standard one of differential geometry, and it is carried out
with [Symbolics.jl](https://symbolics.juliasymbolics.org/):

1. The tangent map ``DF`` is the derivative of the chart; its inverse ``\bar{DF}`` is the
   derivative of the inverse chart, pulled back onto ``\xi``.
2. The metric is inverted symbolically, and both are differentiated twice.
3. The magnetic two-form is the exterior derivative of the vector potential, and the magnetic
   field is its Hodge dual. The Hodge star is orientation-dependent and is given the *signed*
   determinant ``\det DF = \mathrm{orientation} \cdot J``, not the volume element.
4. The magnitude, the unit vector and the perpendicular frame follow, the frame being built from
   the first coordinate basis vector whose cross product with ``b`` does not vanish identically.
5. The electric field is minus the gradient of the potential.
6. Everything is differentiated as far as the tables in [Interface](interface.md) say.

The full derivation is in [Fields](fields.md).

Each resulting expression is then handed to `Symbolics.build_function`, which turns it into
Julia code and compiles it. The functions are `RuntimeGeneratedFunction`s: ordinary compiled
Julia functions that were built at run time rather than parsed from a file. They are stored in
the `FieldFunctions` object, which is the reason a field is a value that can be passed around
rather than a set of names in a module.


## Common Subexpression Elimination

Symbolic differentiation produces expressions that share a great deal of structure and state none
of it. The second derivatives of the Solov'ev flux function contain the same logarithm dozens of
times over, and written out literally they are both unreadable and slow.

`build_function` is therefore asked to eliminate common subexpressions: each repeated
subexpression is computed once into a local and referred to afterwards. This is on by default and
is value-preserving — subexpressions are named, never rewritten — and it can be turned off:

```julia
field = FieldFunctions(equ; cse = false)
```

The only reason to do so is to read the generated code against a paper. The result is
considerably slower.


## What It Costs

The trace, the code generation and the first compilation all happen inside the call to
`FieldFunctions`, and they are not free. A simple cartesian field takes a few seconds; the
Solov'ev equilibrium with an X-point, whose flux function is the largest expression in the
package, takes several. Most of that is Julia compiling the generated code, and it happens once
per field per session.

Two consequences are worth planning around:

* **Build the field once and reuse it.** Constructing it inside a loop rebuilds and recompiles
  everything each time.
* **The cost is per session.** Because the functions are generated at run time, they cannot be
  cached into a precompiled package image. A package that wants a field at load time pays the
  build cost every time it loads.

Evaluating a field, by contrast, is cheap: every accessor is allocation-free and type-stable, and
whole tensors are computed in one pass with their subexpressions shared, which is usually faster
than computing the components one at a time.


## Adding a Field

Adding a field to the package, or defining one in your own code, means adding methods to the
generics in the first table. For a field in cartesian coordinates, the vector potential is
enough:

```@example generation
using ElectromagneticFields
import ElectromagneticFields: A₁, A₂, A₃, get_functions, X, Y, Z

struct MyPinch{T <: Number} <: ElectromagneticFields.CartesianEquilibrium
    name::String
    B₀::T
    MyPinch(B₀::T) where {T} = new{T}("MyPinch", B₀)
end

A₁(x::AbstractVector, equ::MyPinch) = -equ.B₀ * Y(x, equ) / 2
A₂(x::AbstractVector, equ::MyPinch) = +equ.B₀ * X(x, equ) / 2
A₃(x::AbstractVector, equ::MyPinch) = zero(eltype(x))

get_functions(::MyPinch) = (X = X, Y = Y, Z = Z)

field = FieldFunctions(MyPinch(2.0))
B♭(field, 0.0, [0.5, 0.5, 0.5])
```

Three rules govern what those methods may contain.

**Write them generically.** The argument is annotated `AbstractVector` and the element type is
left open, because the generator calls them with a vector of symbolic variables. Annotating
`AbstractVector{Float64}`, or calling a function that only accepts numbers, fails at trace time.
`zero(eltype(x))` rather than `0.0` is the same rule applied to constants.

**Every function used must be differentiable by Symbolics.** The elementary functions are, as are
`NaNMath.log` and the two-argument `atan`, both of which this package relies on. A branch on the
value of a coordinate is not: `if x[1] > 0` cannot be traced, and the expression has to be written
without it.

**A new chart needs its geometry to be consistent.** Defining `x¹` … `ξ³`, the metric, `J` and
`orientation` is not enough on its own; they have to fit together. The test suite asserts
``J = \sqrt{\det(DF^T DF)}`` and ``\det DF = \mathrm{orientation} \cdot J`` for every equilibrium,
and the identities in [Interface](interface.md) are the same checks to run on a chart of your own.
Handedness in particular has no visible symptom when it is wrong other than a magnetic field
pointing the wrong way — see [`orientation`](@ref ElectromagneticFields.orientation).

A perturbation is defined the same way, subtyping `AnalyticPerturbation` (or
`CartesianPerturbation`) and defining whichever of `A₁`, `A₂`, `A₃` and `φ` it contributes. It
must live on the same chart as the equilibrium it perturbs, and the generator asserts that the two
agree on `J` and `g`. The two are added together symbolically before any code is generated, so
`FieldFunctions(equ, pert)` is a field like any other and costs nothing extra to evaluate.


## When It Goes Wrong

**A `MethodError` during `FieldFunctions`** almost always means one of the methods above is not
generic enough — a concrete type annotation, or a call into something that does not accept
symbolic arguments. The stack trace names the method.

**A field that evaluates to `NaN`** on a chart with a singular point is usually the chart rather
than the generator: `1/cos²θ` gauges and `log R` potentials are genuinely undefined where they are
undefined. `NaNMath` variants are used where a domain error would otherwise propagate.

**A magnetic field pointing the wrong way** is an orientation error. Compare against a
finite-difference curl of the vector potential in cartesian coordinates, which is what
`test_curl` in the test suite does.

**A disagreement with an earlier version** can be checked against
`scripts/verify_against_symengine.jl`, which evaluates every quantity for all twenty equilibria
and compares them with values recorded from the SymEngine-based implementation this one replaced.
