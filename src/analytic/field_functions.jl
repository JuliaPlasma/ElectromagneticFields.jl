
"""
    FieldFunction(f)

Wraps a function built by `Symbolics.build_function` and fixes two things about it.

`build_function` builds its output container `similarto` the argument container, so the same
generated function returns an `SMatrix` for an `SVector` argument and a heap `Matrix` for a
`Vector` one — the latter allocating, up to 6 KiB for a rank-three tensor. The wrapper converts the
coordinates to an `SVector` first, so every call is allocation free whatever the caller passes.

It also converts the result to `float(eltype(ξ))`. A structurally constant body — `g♭` of a
cartesian chart, `φ` and `E♭` of a purely magnetic equilibrium — is emitted with `Int` literals and
would otherwise return `Int`, making every component of one field disagree with the others and
forcing a runtime promotion in the caller.
"""
struct FieldFunction{F}
    f::F
end

@inline _astype(::Type{T}, x::Number) where {T} = convert(T, x)
@inline function _astype(::Type{T}, x::StaticArray) where {T}
    convert(similar_type(typeof(x), T), x)
end

@inline function (fun::FieldFunction)(t, ξ::AbstractVector)
    _astype(float(eltype(ξ)), fun.f(t, SVector{3}(ξ)))
end
@inline (fun::FieldFunction)(t, ξ₁, ξ₂, ξ₃) = fun(t, SVector(ξ₁, ξ₂, ξ₃))

"""
The names of the generated functions a [`FieldFunctions`](@ref) holds, which are also the keys of
`functions(field)` and the names of the accessors exported for them.
"""
const FIELD_FUNCTION_NAMES = (
    # chart
    :to_cartesian, :from_cartesian, :DF, :DF̄, :J, :rangemin, :rangemax,
    # metric
    :g♭, :g♯, :Dg♭, :Dg♯, :DDg♭, :DDg♯,
    # vector potential and scalar potential
    :A♭, :A♯, :DA♭, :DDA♭, :φ,
    # magnetic field
    :B, :DB, :DDB, :B♭, :B♯, :B♮, :B♭♭, :DB♭,
    # unit magnetic field
    :b♭, :b♯, :b♮, :Db♭, :Db♮, :DDb♭,
    # perpendicular frame
    :a♭, :a♯, :a♮, :c♭, :c♯, :c♮,
    # electric field
    :E♭, :E♯, :DE♭
)

@doc raw"""
    FieldFunctions(equ, pert = ZeroPerturbation(); cse = true)

The electromagnetic field of `equ`, perturbed by `pert`, as a value: a struct of compiled functions
built from a symbolic trace of the equilibrium's chart, metric and vector potential.

Unlike a set of functions evaluated into a module, this is an ordinary object. It can be built
inside a function, stored, and passed to a vector field that depends on the field — typically in a
`GeometricEquations` problem's `parameters`, from which the equation function evaluates the
components it needs.

# Accessors

Every stored function is reached through an exported generic taking the field first, and accepts
either a coordinate vector or three scalars:

```julia
field = FieldFunctions(ThetaPinchEquilibrium())
B♭(field, t, ξ)            # covariant components of B, an SVector{3}
B♭(field, t, ξ₁, ξ₂, ξ₃)   # the same
B♭(field, t, ξ)[1]         # one component
B(field, t, ξ)             # |B|, a scalar
```

A quantity that has several representations carries the musical isomorphism that produces it:
`♭` for a lowered index, `♯` for a raised one, `♮` for physical components. The bare letter is the
magnitude, so `B` is `|B|` and `DB`, `DDB` are its gradient and Hessian.

| | covariant | contravariant | physical | derivatives |
|:--|:--|:--|:--|:--|
| vector potential | `A♭` | `A♯` | | `DA♭`, `DDA♭` |
| magnetic field | `B♭` | `B♯` | `B♮` | `DB♭`; two-form `B♭♭`; magnitude `B`, `DB`, `DDB` |
| unit magnetic field | `b♭` | `b♯` | `b♮` | `Db♭`, `Db♮`, `DDb♭` |
| perpendicular frame | `a♭`, `c♭` | `a♯`, `c♯` | `a♮`, `c♮` | |
| electric field | `E♭` | `E♯` | | `DE♭` |
| metric | `g♭` | `g♯` | | `Dg♭`, `Dg♯`, `DDg♭`, `DDg♯` |

The chart contributes `to_cartesian`, `from_cartesian`, `DF`, `DF̄`, `J`, `rangemin` and `rangemax`,
and the scalar potential `φ`.

Four things are data rather than functions: `parameters(field)` returns the equilibrium's scalar
parameters, `coordinates(field)` the equilibrium's own coordinate helpers such as `R`, `r`, `θ` and
`ϕ`, `periodicity(field)` the periodic domain, and `orientation(field)` the sign of the chart's
handedness — see [`orientation`](@ref). `functions(field)` returns all of the above generated
functions as a `NamedTuple`.

`cse` names each repeated subexpression once instead of emitting it in full every time it occurs.
It is on by default and is value-preserving; turning it off makes the generated code easier to read
against a paper and considerably slower.
"""
struct FieldFunctions{ET, PT, PAR <: NamedTuple, CRD <: NamedTuple, PER, FNS <: NamedTuple}
    equilibrium::ET
    perturbation::PT

    parameters::PAR
    coordinates::CRD
    periodicity::PER
    orientation::Int

    functions::FNS
end

_outofplace(f::Tuple) = f[1]
_outofplace(f) = f

function FieldFunctions(equ::AnalyticEquilibrium,
        pert::AnalyticPerturbation = ZeroPerturbation(); cse = true)
    exprs = generate_field_expressions(equ, pert)
    t, ξ = exprs.arguments

    function build(ex)
        FieldFunction(_outofplace(Symbolics.build_function(
            ex, t, ξ; expression = Val{false}, cse = cse)))
    end

    fns = NamedTuple{FIELD_FUNCTION_NAMES}(
        map(name -> build(getfield(exprs, name)), FIELD_FUNCTION_NAMES))

    crd = NamedTuple{keys(exprs.coordinates)}(map(build, values(exprs.coordinates)))

    names = if hasmethod(get_parameters, Tuple{typeof(equ)})
        Tuple(get_parameters(equ))
    else
        Tuple(name for name in fieldnames(typeof(equ)) if name != :name)
    end
    par = NamedTuple{names}(map(name -> getfield(equ, name), names))

    FieldFunctions(equ, pert, par, crd,
        GeometricBase.periodicity(zeros(3), equ), orientation(equ), fns)
end

# `FIELD_FUNCTION_NAMES` drives both the construction above and the accessors below, so the two
# cannot drift apart. The one thing that can is the doctable in the docstring.
for name in FIELD_FUNCTION_NAMES
    @eval begin
        @inline $name(field::FieldFunctions, t, ξ::AbstractVector) = field.functions.$name(t, ξ)
        @inline function $name(field::FieldFunctions, t, ξ₁, ξ₂, ξ₃)
            field.functions.$name(t, ξ₁, ξ₂, ξ₃)
        end
    end
end

# The domain bounds depend on neither t nor ξ for every equilibrium in this package, so they are
# also reachable without a time.
@inline rangemin(field::FieldFunctions, ξ::AbstractVector) = rangemin(field, 0, ξ)
@inline rangemax(field::FieldFunctions, ξ::AbstractVector) = rangemax(field, 0, ξ)
@inline rangemin(field::FieldFunctions, ξ₁, ξ₂, ξ₃) = rangemin(field, 0, ξ₁, ξ₂, ξ₃)
@inline rangemax(field::FieldFunctions, ξ₁, ξ₂, ξ₃) = rangemax(field, 0, ξ₁, ξ₂, ξ₃)

"""
    coordinates(field)

The equilibrium's own coordinate helpers as a `NamedTuple` of functions — `X`, `Y`, `Z`, `R`, `r`,
`θ`, `ϕ` and `r²` for the equilibria that define them. Which names are present varies per
equilibrium, which is why these are not accessors of their own.
"""
coordinates(field::FieldFunctions) = field.coordinates

GeometricBase.functions(field::FieldFunctions) = field.functions
GeometricBase.parameters(field::FieldFunctions) = field.parameters
GeometricBase.periodicity(field::FieldFunctions) = field.periodicity
orientation(field::FieldFunctions) = field.orientation

equilibrium(field::FieldFunctions) = field.equilibrium
perturbation(field::FieldFunctions) = field.perturbation

function Base.show(io::IO, field::FieldFunctions)
    print(io, "FieldFunctions for\n")
    print(io, field.equilibrium)
    if !(field.perturbation isa ZeroPerturbation)
        print(io, "\n   perturbed by \n")
        print(io, field.perturbation)
    end
    print(io, "\n Parameters: \n   ", field.parameters)
    print(io, "\n Coordinates: \n   ", keys(field.coordinates))
    print(io, "\n Orientation: ", field.orientation)
end
