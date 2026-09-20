
"""
    FieldFunction(f, p)

A function built by `Symbolics.build_function` together with the parameter values to call it with.

The generated code is free of the equilibrium's parameters: it takes them as an argument, so one
compiled function serves every parameter value of an equilibrium type. `p` holds the values from
the equilibrium struct, flattened by [`parameter_values`](@ref), and the wrapper supplies them so
that callers do not have to.

The wrapper also fixes two things about the generated function. `build_function` builds its output
container `similarto` the argument container, so the same generated function returns an `SMatrix`
for an `SVector` argument and a heap `Matrix` for a `Vector` one — the latter allocating, up to
6 KiB for a rank-three tensor. Converting the coordinates to an `SVector` first keeps every call
allocation free whatever the caller passes.

And it converts the result to `float(eltype(ξ))`. A structurally constant body — `g♭` of a
cartesian chart, `φ` and `E♭` of a purely magnetic equilibrium — is emitted with `Int` literals and
would otherwise return `Int`, making every component of one field disagree with the others and
forcing a runtime promotion in the caller.
"""
struct FieldFunction{F, P}
    f::F
    p::P
end

@inline _astype(::Type{T}, x::Number) where {T} = convert(T, x)
@inline function _astype(::Type{T}, x::StaticArray) where {T}
    convert(similar_type(typeof(x), T), x)
end

@inline function (fun::FieldFunction)(t, ξ::AbstractVector)
    _astype(float(eltype(ξ)), fun.f(t, SVector{3}(ξ), fun.p))
end
@inline (fun::FieldFunction)(t, ξ₁, ξ₂, ξ₃) = fun(t, SVector(ξ₁, ξ₂, ξ₃))

"""
Every quantity a [`FieldFunctions`](@ref) generates, as `name => (shape, description)`.

This is the single source for three things that must not drift apart: the fields of the generated
`NamedTuple`, the accessors exported for them, and their docstrings. `x` in a description is the
cartesian coordinates and `ξ` the chart's own; derivatives are always with respect to `ξ`.
"""
const FIELD_FUNCTION_DOCS = (
    # chart
    to_cartesian = ("`SVector{3}`", "The chart map, ``x(ξ)``."),
    from_cartesian = ("`SVector{3}`", "The inverse chart map, ``ξ(x)``."),
    DF = ("`SMatrix{3,3}`", "The tangent map ``{DF^i}_j = ∂x^i / ∂ξ^j``."),
    DF̄ = ("`SMatrix{3,3}`", "The inverse tangent map ``∂ξ^i / ∂x^j``."),
    J = ("scalar",
        "The volume element ``\\sqrt{|g|} = |\\det DF|``. It carries no sign; see " *
        "[`orientation`](@ref ElectromagneticFields.orientation)."),
    rangemin = ("`SVector{3}`", "Lower bounds of the coordinate domain."),
    rangemax = ("`SVector{3}`", "Upper bounds of the coordinate domain."),

    # metric
    g♭ = ("`SMatrix{3,3}`", "The metric ``g_{ij}``."),
    g♯ = ("`SMatrix{3,3}`", "The inverse metric ``g^{ij}``."),
    Dg♭ = ("`SArray{3,3,3}`", "``∂_k g_{ij}``, indexed `[i,j,k]`."),
    Dg♯ = ("`SArray{3,3,3}`", "``∂_k g^{ij}``, indexed `[i,j,k]`."),
    DDg♭ = ("`SArray{3,3,3,3}`", "``∂_l ∂_k g_{ij}``, indexed `[i,j,k,l]`."),
    DDg♯ = ("`SArray{3,3,3,3}`", "``∂_l ∂_k g^{ij}``, indexed `[i,j,k,l]`."),

    # vector potential and scalar potential
    A♭ = ("`SVector{3}`",
        "Covariant components ``A_i`` of the vector potential, as the " *
        "equilibrium supplied them."),
    A♯ = ("`SVector{3}`", "Contravariant components ``A^i = g^{ij} A_j``."),
    DA♭ = ("`SMatrix{3,3}`", "``∂_j A_i``, indexed `[i,j]`."),
    DDA♭ = ("`SArray{3,3,3}`", "``∂_k ∂_j A_i``, indexed `[i,j,k]`."),
    φ = ("scalar", "The electrostatic potential, zero unless the field defines one."),

    # magnetic field
    B = ("scalar", "The magnitude ``|B| = \\sqrt{B^i B_i}``."),
    DB = ("`SVector{3}`", "The gradient ``∂_i |B|``."),
    DDB = ("`SMatrix{3,3}`", "The Hessian ``∂_j ∂_i |B|``."),
    B♭ = ("`SVector{3}`", "Covariant components ``B_i`` of the magnetic field."),
    B♯ = ("`SVector{3}`", "Contravariant components ``B^i``."),
    B♮ = ("`SVector{3}`", "Physical components, in the orthonormal frame."),
    B♭♭ = ("`SMatrix{3,3}`",
        "The magnetic two-form ``B_{ij} = (∂_i A_j - ∂_j A_i)/2``, " *
        "antisymmetric by construction."),
    DB♭ = ("`SMatrix{3,3}`", "``∂_j B_i``, indexed `[i,j]`."),

    # unit magnetic field
    b♭ = ("`SVector{3}`", "Covariant components of the unit vector ``b = B/|B|``."),
    b♯ = ("`SVector{3}`", "Contravariant components of ``b``."),
    b♮ = ("`SVector{3}`", "Physical components of ``b``."),
    Db♭ = ("`SMatrix{3,3}`", "``∂_j b_i``, indexed `[i,j]`."),
    Db♮ = ("`SMatrix{3,3}`", "``∂_j`` of the physical components of ``b``."),
    DDb♭ = ("`SArray{3,3,3}`", "``∂_k ∂_j b_i``, indexed `[i,j,k]`."),

    # perpendicular frame
    a♭ = (
        "`SVector{3}`", "Covariant components of the first vector perpendicular to ``b``."),
    a♯ = ("`SVector{3}`", "Contravariant components of ``a``."),
    a♮ = ("`SVector{3}`", "Physical components of ``a``."),
    c♭ = ("`SVector{3}`", "Covariant components of ``c = b × a``, completing the triad."),
    c♯ = ("`SVector{3}`", "Contravariant components of ``c``."),
    c♮ = ("`SVector{3}`", "Physical components of ``c``."),

    # electric field
    E♭ = ("`SVector{3}`", "Covariant components ``E_i = -∂_i φ`` of the electric field."),
    E♯ = ("`SVector{3}`", "Contravariant components ``E^i``."),
    DE♭ = ("`SMatrix{3,3}`", "``∂_j E_i``, indexed `[i,j]`.")
)

"""
The names of the generated functions a [`FieldFunctions`](@ref) holds, which are also the keys of
`functions(field)` and the names of the accessors exported for them.
"""
const FIELD_FUNCTION_NAMES = keys(FIELD_FUNCTION_DOCS)

@doc raw"""
    FieldFunctions(equ, pert = ZeroPerturbation(); cse = true, cache = true, cache_module)

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

# Parameters

The generated code is free of the equilibrium's parameters: it is traced in terms of `R₀`, `B₀`,
`q₀` … and takes their values as an argument, which the field supplies from the equilibrium struct
on every call. Two fields of the same type therefore share their compiled code exactly, and
differ only in the values they pass — which is what lets one precompiled function set serve every
parameter value. See `parameter_names` for what counts as a parameter.

`cache` reuses the generated functions of a type already built, which is what makes constructing
any of the shipped equilibria free — they are traced during this package's precompilation and the
cache survives into its image. Redefining a method a field is built from, under Revise say, makes
the cached code stale; `clear_field_cache!` discards it. `cache_module` says where the generated
bodies are kept, and a package building a field during its own precompilation must set it; see
[`@precompilable_fields`](@ref).

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

"""
    @precompilable_fields

Prepare the calling module to hold fields built during its own precompilation.

A generated function keeps its body in a cache inside some module, and that cache only survives
precompilation when the module is the one being precompiled. Put this at the top level of a
package that builds a field there, and pass that module to [`FieldFunctions`](@ref) as
`cache_module`:

```julia
module MyModel

using ElectromagneticFields
@precompilable_fields

const FIELD = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 2.0);
                             cache_module = @__MODULE__)
end
```

`MyModel` then carries the field in its package image: no trace, no code generation and no
compilation happen when it is loaded. Add a `PrecompileTools.@compile_workload` calling the
accessors to cache their compiled code as well, and the first evaluation in a session costs
microseconds.

Without this the field still works, but it is rebuilt on every load.
"""
macro precompilable_fields()
    esc(:($RuntimeGeneratedFunctions.init($__module__)))
end

"""
The generated functions built so far, keyed by what determines them.

Because the parameters are arguments rather than literals, the code depends on the equilibrium's
type and not on its values, so every field of a type already built is a lookup: no trace, no code
generation and no compilation. The entries hold the bare generated functions; the parameter values
are attached per field when it is constructed.

Populated during this package's precompilation for the equilibria in its workload, and it survives
into the package image, so those types cost nothing in a fresh session.

The cache assumes the methods defining a field do not change within a session. Redefine an `A₁` or
a metric coefficient — under Revise, say — and [`clear_field_cache!`](@ref) is what makes the next
`FieldFunctions` see it.
"""
const FIELD_CACHE = Dict{Any, Any}()

const FIELD_CACHE_LOCK = ReentrantLock()

"""
    clear_field_cache!()

Discard the generated functions cached by [`FieldFunctions`](@ref), so that the next field of each
type is traced afresh. Needed only after redefining a method that a field is built from.
"""
function clear_field_cache!()
    @lock FIELD_CACHE_LOCK empty!(FIELD_CACHE)
    nothing
end

# The trace and the code generation, without the parameter values.
function _generate_functions(equ, pert, cse, cache_module)
    exprs = generate_field_expressions(equ, pert)
    t, ξ, p = exprs.arguments

    build(ex) = _outofplace(Symbolics.build_function(
        ex, t, ξ, p; expression = Val{false}, cse = cse,
        expression_module = cache_module))

    (
        functions = NamedTuple{FIELD_FUNCTION_NAMES}(
            map(name -> build(getfield(exprs, name)), FIELD_FUNCTION_NAMES)),
        coordinates = map(build, exprs.coordinates),
        nparameters = length(p))
end

function FieldFunctions(equ::AnalyticEquilibrium,
        pert::AnalyticPerturbation = ZeroPerturbation();
        cse = true, cache_module = @__MODULE__, cache = true)
    # The values the generated code is called with, taken from the equilibrium and perturbation
    # structs and flattened in the same order as the symbols they replace. An `SVector` so that
    # passing them costs nothing.
    pvalues = SVector(parameter_values(equ)..., parameter_values(pert)...)

    # What the generated code depends on, and nothing more. The parameter count is in the key
    # because a parameter may itself be a vector — the Solov'ev coefficients — whose length a type
    # does not by itself pin down.
    key = (ConstructionBase.constructorof(typeof(equ)),
        ConstructionBase.constructorof(typeof(pert)),
        length(pvalues), cse, cache_module)

    generated = if cache
        @lock FIELD_CACHE_LOCK get!(
            () -> _generate_functions(equ, pert, cse, cache_module), FIELD_CACHE, key)
    else
        _generate_functions(equ, pert, cse, cache_module)
    end

    @assert generated.nparameters == length(pvalues)

    fns = map(f -> FieldFunction(f, pvalues), generated.functions)
    crd = map(f -> FieldFunction(f, pvalues), generated.coordinates)

    names = parameter_names(equ)
    par = NamedTuple{names}(map(name -> getfield(equ, name), names))

    FieldFunctions(equ, pert, par, crd,
        GeometricBase.periodicity(zeros(3), equ), orientation(equ), fns)
end

# `FIELD_FUNCTION_DOCS` drives the construction above, the accessors below and their docstrings,
# so the three cannot drift apart.
for name in FIELD_FUNCTION_NAMES
    shape, description = FIELD_FUNCTION_DOCS[name]

    docstring = """
        $name(field::FieldFunctions, t, ξ)
        $name(field::FieldFunctions, t, ξ₁, ξ₂, ξ₃)

    $description

    Returns $shape, of `float(eltype(ξ))`. Allocation free for any coordinate container.

    See [`FieldFunctions`](@ref) for the naming scheme and the full list of quantities.
    """

    @eval begin
        @doc $docstring @inline function $name(
                field::FieldFunctions, t, ξ::AbstractVector)
            field.functions.$name(t, ξ)
        end

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
