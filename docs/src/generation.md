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
| `get_parameters` | which fields are parameters of the field | every field but `name` |

A field built on `CartesianEquilibrium` inherits the identity chart, `J = 1`, the euclidean
metric and `orientation = +1`, so in that case only the vector potential is left. That is why
`ThetaPinchEquilibrium` is a few dozen lines.

The struct itself holds the parameters — `B₀` for the θ-pinch, `R₀`, `B₀` and `q₀` for the
tokamaks — and the methods read them off it. They are not baked into the generated code; see
[The Parameters Are Arguments, Not Literals](@ref).


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


## The Parameters Are Arguments, Not Literals

One detail of the trace decides a great deal. It does not run against the equilibrium you passed
in, but against a copy of it whose parameters are symbolic. `AxisymmetricTokamakCylindrical`
is traced in terms of ``R_0``, ``B_0`` and ``q_0`` rather than `6.2`, `5.3` and `2.0`, so the
generated code takes them as an argument. The equilibrium struct remains where the values live,
and the field passes them in on every call.

The consequence is that **the generated code depends on the equilibrium's type, not on its
parameters**. Two fields of the same type share it exactly:

```@example generation
using ElectromagneticFields

a = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(1.0, 1.0, 2.0))
b = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 1.7))

shared = typeof(functions(a).B♭.f) === typeof(functions(b).B♭.f)
@assert shared # hide
shared, B♭(a, 0.0, [1.05, 0.25, 0.5]), B♭(b, 0.0, [1.05, 0.25, 0.5])
```

Same code, different answers. This is why the precompilation described below works for parameter
values nobody anticipated, and it is why `parameter_names` matters: anything the `A₁`, `φ` or
metric methods read that is *not* listed as a parameter gets frozen into the code as a literal.
The default — every field of the struct except `name` — is what the equilibria here want.


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


## Reading the Generated Code

There is no equivalent of printing the generated module, because there is no module. To see what
was generated for one quantity, take the symbolic expressions and ask `build_function` for an
expression instead of a function:

```@example generation
using ElectromagneticFields: generate_field_expressions
using Symbolics

exprs = generate_field_expressions(AxisymmetricTokamakCylindricalEquilibrium(),
                                   ZeroPerturbation())
t, ξ, p = exprs.arguments

Symbolics.build_function(exprs.B♭, t, ξ, p; expression = Val{true}, cse = true)
```

`exprs.arguments` are the symbols the code is written in: the time, the three coordinates, and
the parameters — which appear in the result as `R₀`, `B₀` and `q₀` rather than as numbers, since
the code takes them as its third argument. The rest of `exprs` is a `NamedTuple` keyed exactly
like `functions(field)`, so any quantity in the tables of [Interface](interface.md) can be
inspected this way. Passing `cse = false` gives the unshared form, which is the one to read
against a paper.


## What It Costs

The trace, the code generation and the first compilation all happen inside the call to
`FieldFunctions`. Evaluating the result afterwards is cheap — every accessor is allocation-free
and type-stable, and a whole tensor is computed in one pass with its subexpressions shared, which
is usually faster than computing the components one at a time.

Building it is the part worth understanding. Measured on the toroidal tokamak in a cold session:

| stage | time |
|:--|--:|
| symbolic trace | 7.1 s |
| code generation (≈8400 lines of Julia) | 4.1 s |
| constructing the generated functions | 0.2 s |
| compiling the generated code (first call) | 0.6 s |

The striking thing is where the time is *not*. Compiling the generated code is 5% of the bill.
Almost all of it is Julia compiling **Symbolics' own machinery** for the expression types a trace
produces — and that is paid per *expression shape*, not per field. Because the parameters are
arguments, the shape is a property of the equilibrium's type, so building the same equilibrium
again with different parameters costs nothing measurable.

That is what makes the cost recoverable, and this package recovers all of it for the equilibria
it ships. Two things combine:

* the generated functions of a type already built are **cached**, keyed by what actually
  determines them — the equilibrium and perturbation types, the parameter count, `cse` and
  `cache_module` — but not by the parameter values, which the code does not contain;
* a `PrecompileTools` workload **traces one field of every shipped type** during precompilation,
  so both the compiled specializations and the cache itself land in the package image.

A fresh session therefore finds them already built:

| | before | after |
|:--|--:|--:|
| first field built in a session | 8.5 s | 0.00 s |
| all twenty shipped equilibria | ~19 s | 0.01 s |
| any parameter value of those types | — | 0.00 s |

at the price of this package's own precompilation, 1.5 s → 24.5 s, paid once per version.

An equilibrium you define yourself is traced the first time and cached thereafter, so it costs a
fraction of a second once per session — or nothing at all, if you precompile it in your own
package as described below.

One consequence of the cache is worth knowing. It is a `Dict`, so `FieldFunctions(equ)` is not
type-inferable at its call site: the object it returns is concretely typed and every accessor on
it is type-stable and allocation-free, but the construction itself is not. Fields are normally
built in a constructor or a setup step and then handed to something else, which is a function
barrier already, so this costs nothing in practice — but build a field and evaluate it in the
*same* function body and the evaluation will be inferred as `Any`.

!!! warning "The cache assumes your field definitions do not change"
    It is keyed on types, not on the content of the methods. If you redefine an `A₁` or a metric
    coefficient in a running session — under Revise, say — the next `FieldFunctions` of that type
    returns the code built from the *old* definition. Call `clear_field_cache!()` after such an
    edit. `FieldFunctions(equ; cache = false)` bypasses the cache for one call.


## Precompiling a Field in Your Own Package

The remaining fraction of a second per field can be removed too. A package that always uses the
same field can build it during its *own* precompilation, so loading it costs nothing at all.

Two things are needed. A generated function keeps its body in a cache inside some module, and
that cache only survives precompilation when the module is the one being precompiled;
[`@precompilable_fields`](@ref) prepares the calling module, and `cache_module` tells
`FieldFunctions` to use it. A `PrecompileTools.@compile_workload` then caches the compiled code
of the accessors as well.

```julia
module MyModel

using ElectromagneticFields
using PrecompileTools

@precompilable_fields

const FIELD = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 2.0);
                             cache_module = @__MODULE__)

sample(t, ξ) = B♭(FIELD, t, ξ)

@setup_workload begin
    ξ = [6.5, 0.5, 0.25]
    @compile_workload begin
        sample(0.0, ξ)
    end
end

end
```

Loading `MyModel` in a fresh session then costs only what loading ElectromagneticFields costs, and
the first evaluation of the field takes **microseconds**: no trace, no code generation, no
compilation. Omitting `cache_module` is not an error — the field is simply rebuilt on every load.


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
    MyPinch{T}(B₀::T) where {T <: Number} = new("MyPinch", B₀)
end

MyPinch(B₀::T) where {T <: Number} = MyPinch{T}(B₀)

A₁(x::AbstractVector, equ::MyPinch) = -equ.B₀ * Y(x, equ) / 2
A₂(x::AbstractVector, equ::MyPinch) = +equ.B₀ * X(x, equ) / 2
A₃(x::AbstractVector, equ::MyPinch) = zero(eltype(x))

get_functions(::MyPinch) = (X = X, Y = Y, Z = Z)

field = FieldFunctions(MyPinch(2.0))
B♭(field, 0.0, [0.5, 0.5, 0.5])
```

Three rules govern what those methods may contain.

**Follow the constructor convention.** The inner constructor takes the parameters in order and
sets `name` itself, and there is an outer one supplying defaults — which is what lets the trace
build the copy with symbolic parameters. `MyPinch{T}(B₀::T)`, not `MyPinch(B₀::T)`: the trace
calls `MyPinch{Num}(...)`. A type that cannot be written this way needs a `symbolic_copy` method.

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
