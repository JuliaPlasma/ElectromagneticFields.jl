
using Combinatorics
using LinearAlgebra
using StaticArrays
using Symbolics

abstract type AnalyticField <: ElectromagneticField end
abstract type AnalyticEquilibrium <: AnalyticField end
abstract type AnalyticPerturbation <: AnalyticField end

"""
    get_functions(equ)

The equilibrium's own coordinate helpers, as a `NamedTuple` mapping each name to a method taking
`(x, equ)` — `(X = X, Y = Y, Z = Z, R = R, r = r, θ = θ, ϕ = ϕ, r² = r²)` for the axisymmetric
equilibria. They are generated alongside the field and reached through `coordinates(field)`.

Which names an equilibrium offers is its own business, which is why these are a `NamedTuple`
rather than accessors of their own. Defining no method means the field has none.
"""
function get_functions end

"""
    get_parameters(equ)

The fields of `equ` that are parameters of the electromagnetic field, in the order its
constructor takes them.

Defining no method means every field but `name`, which is what all the equilibria here want.
Define one for a type with a field that is not a parameter — a cache, say — bearing in mind that
the generated code takes exactly these as its argument, so anything the `A₁`, `φ` or metric
methods read and this does not list is frozen into the code as a literal. See
[`parameter_names`](@ref), which applies the default.

A value frozen in that way is invisible to the cache, which keys on the types and the parameter
shapes alone: two equilibria differing only in such a value share one entry, and the second gets
the first one's code. List the value as a parameter, or build the field with `cache = false`.
"""
function get_parameters end

# Coordinate helpers an equilibrium may define for its own chart, each as a method taking
# `(x, equ)`, so one name carries one method per equilibrium. `get_functions` lists which of them
# to generate into `coordinates(field)`. They are deliberately not exported — `r`, `θ` and `ϕ`
# would collide with practically any caller.
function X end
function Y end
function Z end
function R end
function r end
function θ end
function ϕ end
function r² end

function x¹(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("x¹() not implemented for ", ET)
end
function x²(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("x²() not implemented for ", ET)
end
function x³(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("x³() not implemented for ", ET)
end

function ξ¹(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("ξ¹() not implemented for ", ET)
end
function ξ²(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("ξ²() not implemented for ", ET)
end
function ξ³(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("ξ³() not implemented for ", ET)
end

function J(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("J() not implemented for ", ET)
end

@doc raw"""
    orientation(equ)

Sign of the chart's orientation: `+1` if `(ξ¹, ξ², ξ³)` is right-handed, `-1` if left-handed.

`J(x, equ)` is the *volume element* `√|g| = |det DF|`, which is what most of the machinery here
wants and what [`J`](@ref) returns. The Hodge star and the cross product, however, are
orientation-dependent and need the *signed* determinant `det DF = orientation(equ) * J(x, equ)`.

Several of the charts in this package are left-handed — `(R, Z, ϕ)` and `(r, θ, ϕ)` both have
`det DF < 0`, since the right-handed orderings would be `(R, ϕ, Z)` and `(r, ϕ, θ)` — so this is
not a corner case. Getting it wrong reverses `B` without any other visible symptom: the same
vector potential yields a magnetic field antiparallel to the one the cartesian chart gives at the
same physical point.

# Orientation of the charts in this package

Every equilibrium built on `CartesianEquilibrium` uses the identity map `(x, y, z)`, so `DF = I`,
`J = 1` and the chart is right-handed. The four families below define their own coordinates, and
all four are left-handed — in each case because the toroidal angle `ϕ` sits in the third slot where
the right-handed ordering would put the second poloidal coordinate.

| chart | coordinates `(ξ¹, ξ², ξ³)` | `J` | `det DF` | orientation |
|:--|:--|:--|:--|:--:|
| `CartesianEquilibrium` and all its subtypes | `(x, y, z)` | `1` | `+1` | `+1` |
| `AxisymmetricTokamakCylindricalEquilibrium` | `(R, Z, ϕ)` | `R` | `-R` | `-1` |
| `AxisymmetricTokamakToroidalEquilibrium` | `(r, θ, ϕ)` | `r R` | `-r R` | `-1` |
| `AxisymmetricTokamakToroidalRegularizationEquilibrium` | `(r, θ, ϕ)` | `r R` | `-r R` | `-1` |
| `AbstractSolovevEquilibrium` | `(R/R₀, Z/R₀, ϕ)` | `R R₀²` | `-R R₀²` | `-1` |

Concretely: `AxisymmetricTokamakCartesianEquilibrium`, `SolovevSymmetricEquilibrium`,
`ThetaPinchEquilibrium`, `DipoleField`, `ABCEquilibrium`, `SingularEquilibrium`,
`SymmetricQuadraticEquilibrium`, `QuadraticPotentialsField` and the three Penning traps are
right-handed; `AxisymmetricTokamakCylindricalEquilibrium`, `AxisymmetricTokamakToroidalEquilibrium`,
`AxisymmetricTokamakToroidalRegularizationEquilibrium` and every Solov'ev equilibrium other than
`SolovevSymmetricEquilibrium` (which is a cartesian chart despite the name) are left-handed.

[`FieldFunctions`](@ref) stores the sign as a value rather than generating a function for it, since
it depends on neither `t` nor `ξ`, and `orientation(field)` returns it. A loaded field therefore
recovers `det DF = orientation(field) * J(field, t, ξ)` without the equilibrium object.

`test_analytic.jl` asserts `det(DF) ≈ orientation * J` for each of the charts above, and separately
that the stored value equals the trait it was generated from. A new chart that declares the wrong
sign therefore fails immediately rather than silently flipping its own `B`, and so does a generator
that stops tracking the trait.
"""
orientation(::AnalyticField) = 1

g₁₁(::AbstractVector{T}, ::AnalyticField) where {T} = one(T)
g₁₂(::AbstractVector{T}, ::AnalyticField) where {T} = zero(T)
g₁₃(::AbstractVector{T}, ::AnalyticField) where {T} = zero(T)
g₂₁(::AbstractVector{T}, ::AnalyticField) where {T} = zero(T)
g₂₂(::AbstractVector{T}, ::AnalyticField) where {T} = one(T)
g₂₃(::AbstractVector{T}, ::AnalyticField) where {T} = zero(T)
g₃₁(::AbstractVector{T}, ::AnalyticField) where {T} = zero(T)
g₃₂(::AbstractVector{T}, ::AnalyticField) where {T} = zero(T)
g₃₃(::AbstractVector{T}, ::AnalyticField) where {T} = one(T)

minx¹(ξ::AbstractVector{T}, equ::AnalyticField) where {T} = -T(Inf)
minx²(ξ::AbstractVector{T}, equ::AnalyticField) where {T} = -T(Inf)
minx³(ξ::AbstractVector{T}, equ::AnalyticField) where {T} = -T(Inf)

maxx¹(ξ::AbstractVector{T}, equ::AnalyticField) where {T} = +T(Inf)
maxx²(ξ::AbstractVector{T}, equ::AnalyticField) where {T} = +T(Inf)
maxx³(ξ::AbstractVector{T}, equ::AnalyticField) where {T} = +T(Inf)

# `periodic(equ)` says which of the chart's three coordinates are periodic, one `Bool` each. It is
# defined per chart family, in that chart's own source file next to its `minx`/`maxx` bounds, and
# there is deliberately **no default**: a chart nobody has answered for raises a `MethodError` when
# a field is built from it.
#
# The name is GeometricBase's, where `periodic(s::StateVariable)` already means one `Bool` per
# component and `isperiodic` means `any` of them. GeometricBase's other generic, `periodicity`, is
# not this: GeometricEquations gives it an `(xmin, xmax)` tuple, and derives the per-component
# answer from those bounds under the name `getperiodicity`. Answering `periodicity` with a `Bool`
# vector would give one generic two shapes, and `per_lo, per_hi = periodicity(equ)` on a `Bool`
# vector destructures to `(false, false)` without error.
#
# It cannot be derived from `minx`/`maxx`, because a bounded range does not imply periodicity. The
# two agree on every chart here — each bounds exactly its angles, to [0, 2π] — but a wall at r = a
# or a slab bounded in z would bound a coordinate that does not wrap. An all-`false` default fails
# for the same reason a derived answer does: it lets a periodic chart answer silently.

"""
    from_cartesian(x, equ::AnalyticField)

Evaluate the inverse chart map ``ξ(x)`` of an equilibrium directly, without generating any code.

This two-argument method shares its name with the [`FieldFunctions`](@ref) accessor
`from_cartesian(field, t, x)`, and the number of arguments is what picks between them.
`from_cartesian(x, equ)` reaches this method, which returns an allocating `Vector`; the accessor
is the three-argument one and returns an `SVector{3}`. Prefer the accessor unless the equilibrium
is all you have.
"""
from_cartesian(x::AbstractVector, equ::AnalyticField) = [ξ¹(x, equ), ξ²(x, equ), ξ³(x, equ)]

"""
    to_cartesian(ξ, equ::AnalyticField)

Evaluate the chart map ``x(ξ)`` of an equilibrium directly, without generating any code.

This two-argument method shares its name with the [`FieldFunctions`](@ref) accessor
`to_cartesian(field, t, ξ)`, and the number of arguments is what picks between them.
`to_cartesian(ξ, equ)` reaches this method, which returns an allocating `Vector`; the accessor is
the three-argument one and returns an `SVector{3}`. Prefer the accessor unless the equilibrium is
all you have.
"""
to_cartesian(ξ::AbstractVector, equ::AnalyticField) = [x¹(ξ, equ), x²(ξ, equ), x³(ξ, equ)]

function A₁(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("A₁() not implemented for ", ET)
end
function A₂(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("A₂() not implemented for ", ET)
end
function A₃(::AbstractVector, ::ET) where {ET <: AnalyticField}
    error("A₃() not implemented for ", ET)
end

φ(::AbstractVector{T}, ::AnalyticField) where {T} = zero(T)

function A(x, equ)
    [A₁(x, equ), A₂(x, equ), A₃(x, equ)]
end

function g(x, equ)
    [g₁₁(x, equ) g₁₂(x, equ) g₁₃(x, equ);
     g₂₁(x, equ) g₂₂(x, equ) g₂₃(x, equ);
     g₃₁(x, equ) g₃₂(x, equ) g₃₃(x, equ)]
end

struct ZeroPerturbation <: AnalyticPerturbation
    name::String
    ZeroPerturbation() = new("ZeroPerturbation")
end

A₁(::AbstractVector{T}, ::AnalyticPerturbation) where {T} = zero(T)
A₂(::AbstractVector{T}, ::AnalyticPerturbation) where {T} = zero(T)
A₃(::AbstractVector{T}, ::AnalyticPerturbation) where {T} = zero(T)

"Returns the i-th component of the vector corresponding to the one-form α"
function covariant_to_contravariant(α, g̅, i)
    g̅[i, 1] * α[1] + g̅[i, 2] * α[2] + g̅[i, 3] * α[3]
end

"Returns the i-th component of the one-form corresponding to the vector v"
function contravariant_to_covariant(v, g, i)
    g[i, 1] * v[1] + g[i, 2] * v[2] + g[i, 3] * v[3]
end

"Returns the i-th component of the physical coordinate representation of the one-form α"
function covariant_to_physical(α, DF̄, i)
    DF̄[1, i] * α[1] + DF̄[2, i] * α[2] + DF̄[3, i] * α[3]
end

"Returns the i-th component of the physical coordinate representation of the vector v"
function contravariant_to_physical(v, DF, i)
    DF[i, 1] * v[1] + DF[i, 2] * v[2] + DF[i, 3] * v[3]
end

"Returns the m-th component of the one-form corresponding to the two-form β"
function hodge²¹(β, g̅, J, m)
    α = 0

    for i in 1:3
        for j in 1:3
            for k in 1:3
                for l in 1:3
                    α += β[i, j] * g̅[i, k] * g̅[j, l] * levicivita([k, l, m])
                end
            end
        end
    end

    return J * α
end

"Returns the m-th component of the cross-product between the vectors v and w"
function crossproduct(v, w, g̅, J, l)
    u = zero(J)

    for i in 1:3
        for j in 1:3
            for k in 1:3
                u += v[i] * w[j] * g̅[k, l] * levicivita([i, j, k])
            end
        end
    end

    return J * u
end

"Returns the length of the vector v"
function magnitude(v, g)
    l = zero(eltype(v))

    for i in 1:3
        for j in 1:3
            l += v[i] * g[i, j] * v[j]
        end
    end

    return sqrt(l)
end

"Normalises the vector v in the metric g"
function normalize(v, g)
    return v ./ magnitude(v, g)
end

"Normalises the vector v in the metric g"
function normalize!(v, g)
    v ./= magnitude(v, g)
end

"""
    parameter_names(field)

The fields of an equilibrium or perturbation that are parameters of the electromagnetic field, in
the order its constructor takes them.

Every field of the struct except `name` by default, which is what the equilibria here want.
Override [`get_parameters`](@ref) for a type with a field that is not a parameter — a cache, say —
bearing in mind that anything the `A₁`/`φ`/metric methods read must be listed, or the generated
code will have it baked in as a literal rather than taking it as an argument.
"""
function parameter_names(field::AnalyticField)
    if hasmethod(get_parameters, Tuple{typeof(field)})
        Tuple(get_parameters(field))
    else
        Tuple(name for name in fieldnames(typeof(field)) if name != :name)
    end
end

"""
    parameter_values(field)

The parameters of `field` flattened into a `Tuple` of scalars, which is the form the generated
code takes them in. A parameter that is itself a vector — the Solov'ev coefficients — contributes
each of its entries.
"""
function parameter_values(field::AnalyticField)
    vals = Any[]
    for name in parameter_names(field)
        value = getfield(field, name)
        value isa Number ? push!(vals, value) : append!(vals, value)
    end
    Tuple(vals)
end

parameter_values(::ZeroPerturbation) = ()

"""
    parameter_shape(field)

How the parameters of `field` spread over the flattened tuple [`parameter_values`](@ref) returns:
one entry per name, `-1` for a scalar and the length for a vector.

The generated code reads its parameter argument positionally, so which slot carries which meaning
follows from this shape and not from the number of slots — two splits totalling the same, such as
`(2, 3)` and `(3, 2)`, map the slots differently. It is part of the key
[`FieldFunctions`](@ref) looks its cache up by.

A scalar takes `-1` rather than `0` so that it cannot read as a vector of length zero. The shape
then determines the number of slots, which makes it strictly finer than that number: two fields
this tells apart are never served each other's code.
"""
function parameter_shape(field::AnalyticField)
    map(parameter_names(field)) do name
        value = getfield(field, name)
        value isa Number ? -1 : length(value)
    end
end

"""
    symbolic_copy(field)

A copy of `field` whose parameters are symbolic, so that tracing through it produces expressions
in the parameters rather than in their values. Returns the copy and the symbols, flattened in the
order [`parameter_values`](@ref) uses.

The reconstruction assumes the convention every field here follows: a parametric struct whose
first member is `name` and whose inner constructor takes the remaining members in order. A type
that departs from it should add a method.
"""
function symbolic_copy(field::AnalyticField, prefix::Symbol)
    names = parameter_names(field)
    isempty(names) && return field, Num[]

    flat = Num[]
    members = Any[]

    for name in names
        value = getfield(field, name)
        if value isa Number
            s = Symbolics.variable(Symbol(prefix, :_, name))
            push!(flat, s)
            push!(members, s)
        else
            ss = [Symbolics.variable(Symbol(prefix, :_, name, :_, i))
                  for i in eachindex(value)]
            append!(flat, ss)
            push!(members, ss)
        end
    end

    wrapper = ConstructionBase.constructorof(typeof(field))
    constructor = try
        wrapper{Num}
    catch
        error(
            "cannot build a symbolic copy of $(typeof(field)): the trace needs to construct ",
            "it with symbolic parameters, which assumes a struct with a single type parameter ",
            "and an inner constructor taking the parameters in order, as in ",
            "`ThetaPinchEquilibrium{T}(B₀::T)`. Add a `symbolic_copy` method for this type.")
    end

    constructor(members...), flat
end

symbolic_copy(field::ZeroPerturbation, ::Symbol) = (field, Num[])

"""
    generate_field_expressions(equ, pert)

Build the symbolic expression for every quantity [`FieldFunctions`](@ref) stores, returned as a
`NamedTuple` whose keys are that struct's field names.

Each value is a `Num`, or a `StaticArray` of `Num` of the quantity's tensor rank: `SVector{3}` for
a one-form or a vector, `SMatrix{3,3}` for a Jacobian or the metric, `SArray{Tuple{3,3,3}}` and
`SArray{Tuple{3,3,3,3}}` for the higher derivatives. The container matters — `build_function`
builds its output `similarto` what it is given.
"""
function generate_field_expressions(
        equilibrium::AnalyticEquilibrium, perturbation::AnalyticPerturbation)
    # Symbols for time t and coordinates x = (x₁, x₂, x₃), ξ = (ξ₁, ξ₂, ξ₃).
    Symbolics.@variables t x₁ x₂ x₃ ξ₁ ξ₂ ξ₃
    x = [x₁, x₂, x₃]
    ξ = [ξ₁, ξ₂, ξ₃]

    # The trace runs against copies whose parameters are symbolic, so the expressions below — and
    # the code built from them — are in terms of R₀, B₀, q₀ … rather than their values. The
    # values travel separately, and one compiled set therefore serves every parameter value of an
    # equilibrium type.
    equ, equ_params = symbolic_copy(equilibrium, :equ)
    pert, pert_params = symbolic_copy(perturbation, :pert)
    p = vcat(equ_params, pert_params)

    D(f, v) = Symbolics.derivative(f, v)

    # check for compatible metric
    if typeof(pert) != ZeroPerturbation
        @assert isequal(J(x, equ), J(x, pert))
        @assert all(isequal.(g(x, equ), g(x, pert)))
    end

    # cartesian coordinates
    x̂ = [x¹(ξ, equ), x²(ξ, equ), x³(ξ, equ)]

    # curvilinear coordinates
    ξ̂ = [ξ¹(x, equ), ξ²(x, equ), ξ³(x, equ)]

    # ranges
    minx̂ = [minx¹(ξ, equ), minx²(ξ, equ), minx³(ξ, equ)]
    maxx̂ = [maxx¹(ξ, equ), maxx²(ξ, equ), maxx³(ξ, equ)]

    # Jacobian
    DF = [D(x̂[i], ξ[j]) for i in 1:3, j in 1:3]

    DF̄ = [Symbolics.substitute(D(ξ̂[i], x[j]),
              Dict(x₁ => x¹(ξ, equ), x₂ => x²(ξ, equ), x₃ => x³(ξ, equ)))
          for i in 1:3, j in 1:3]

    # obtain metric and invert it
    gmat = g(ξ, equ)
    ginv = inv(gmat)

    # derivatives of metric coefficients
    Dg = [D(gmat[i, j], ξ[k]) for i in 1:3, j in 1:3, k in 1:3]
    Dḡ = [D(ginv[i, j], ξ[k]) for i in 1:3, j in 1:3, k in 1:3]
    DDg = [D(D(gmat[i, j], ξ[k]), ξ[l]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]
    DDḡ = [D(D(ginv[i, j], ξ[k]), ξ[l]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3]

    # volume element
    Jdet = J(ξ, equ)

    # Signed Jacobian determinant det(DF), for the orientation-dependent operations below. See
    # `orientation`.
    Jsgn = orientation(equ) * Jdet

    # obtain vector potential
    A¹ = A(ξ, equ) .+ A(ξ, pert)

    # compute vector potential in contravariant coordinates
    Avec = [covariant_to_contravariant(A¹, ginv, i) for i in 1:3]

    # compute Jacobian and second derivative of vector potential A
    DA = [D(A¹[i], ξ[j]) for i in 1:3, j in 1:3]
    DDA = [D(DA[i, j], ξ[k]) for i in 1:3, j in 1:3, k in 1:3]

    # compute components of magnetic field B
    Bᶜ = [DA[3, 2] - DA[2, 3],
        DA[1, 3] - DA[3, 1],
        DA[2, 1] - DA[1, 2]]

    # compute magnetic field two-form B²
    B² = [0 +Bᶜ[3] -Bᶜ[2];
          -Bᶜ[3] 0 +Bᶜ[1];
          +Bᶜ[2] -Bᶜ[1] 0] .* 1 // 2

    # compute magnetic field one-form B¹ = ⋆B²
    B¹ = [hodge²¹(B², ginv, Jsgn, i) for i in 1:3]

    # compute magnetic field in physical and contravariant coordinates
    Bphys = [covariant_to_physical(B¹, DF̄, i) for i in 1:3]
    Bvec = [covariant_to_contravariant(B¹, ginv, i) for i in 1:3]

    # compute absolute value |B| of B
    Babs = sqrt(transpose(Bvec) * B¹)

    # compute magnetic unit one-form, and its physical and contravariant representations
    b¹ = [B¹[i] / Babs for i in 1:3]
    bphys = [Bphys[i] / Babs for i in 1:3]
    bvec = [Bvec[i] / Babs for i in 1:3]

    # compute Jacobians of B and of the magnetic unit one-form b
    DB = [D(B¹[i], ξ[j]) for i in 1:3, j in 1:3]
    Db = [D(b¹[i], ξ[j]) for i in 1:3, j in 1:3]
    Dbphys = [D(bphys[i], ξ[j]) for i in 1:3, j in 1:3]
    DDb = [D(D(b¹[i], ξ[j]), ξ[k]) for i in 1:3, j in 1:3, k in 1:3]

    # compute first and second derivatives of absolute value of magnetic field
    DBabs = [D(Babs, ξ[j]) for j in 1:3]
    DDBabs = [D(D(Babs, ξ[i]), ξ[j]) for i in 1:3, j in 1:3]

    # compute unit vectors perpendicular to magnetic field
    avec = [Num(0), Num(0), Num(0)]
    for tvec in ([Num(1), Num(0), Num(0)],
        [Num(0), Num(1), Num(0)],
        [Num(0), Num(0), Num(1)])
        avec .= [crossproduct(tvec, bvec, ginv, Jsgn, i) for i in 1:3]
        # `iszero` alone decides whether the cross product vanishes. `simplify` must not be
        # used here: it cancels fractions through a polynomial gcd over `Rational{Int64}` that
        # overflows on the coefficients some of these fields carry, and it leaves state on the
        # shared subexpressions that makes a later rebuild of the same field generate a
        # differently ordered — and so not bitwise equal — function.
        if !all(iszero, avec)
            break
        end
    end
    cvec = [crossproduct(bvec, avec, ginv, Jsgn, i) for i in 1:3]

    normalize!(avec, gmat)
    normalize!(cvec, gmat)

    # compute components of magnetic unit vectors in physical coordinates
    aphys = DF * avec
    cphys = DF * cvec

    # compute components of magnetic unit vectors in covariant coordinates
    a¹ = gmat * avec
    c¹ = gmat * cvec

    # obtain scalar potential
    φ⁰ = φ(ξ, equ) .+ φ(ξ, pert)

    # compute components of electric field E, its Jacobian and its contravariant representation
    E¹ = [D(-φ⁰, ξ[i]) for i in 1:3]
    DE = [D(E¹[i], ξ[j]) for i in 1:3, j in 1:3]
    Evec = [covariant_to_contravariant(E¹, ginv, i) for i in 1:3]

    # the per-equilibrium coordinate helpers, e.g. (X, Y, Z, R, r, θ, ϕ, r²)
    coordinates = if hasmethod(get_functions, Tuple{typeof(equ)})
        fs = get_functions(equ)
        NamedTuple{keys(fs)}(map(f -> f(ξ, equ), values(fs)))
    else
        NamedTuple()
    end

    # `subs` on the inverse chart is needed because ξ¹, ξ², ξ³ are written in terms of x
    tocurvilinear = [Symbolics.substitute(ξ̂[i], Dict(x₁ => ξ₁, x₂ => ξ₂, x₃ => ξ₃))
                     for i in 1:3]

    # `build_function` builds its output `similarto` the container it is given, so every expression
    # goes in as a StaticArray of the right rank.
    sv(v) = SVector{3}(Num.(v))
    sm(m) = SMatrix{3, 3}(Num.(m))
    s3(a) = SArray{Tuple{3, 3, 3}}(Num.(a))
    s4(a) = SArray{Tuple{3, 3, 3, 3}}(Num.(a))

    (
        # symbols the caller needs to build functions of these expressions
        arguments = (t, ξ, p),
        coordinates = coordinates,

        # chart
        to_cartesian = sv(x̂),
        from_cartesian = sv(tocurvilinear),
        DF = sm(DF),
        DF̄ = sm(DF̄),
        J = Num(Jdet),
        rangemin = sv(minx̂),
        rangemax = sv(maxx̂),

        # metric
        g♭ = sm(gmat),
        g♯ = sm(ginv),
        Dg♭ = s3(Dg),
        Dg♯ = s3(Dḡ),
        DDg♭ = s4(DDg),
        DDg♯ = s4(DDḡ),

        # vector potential and scalar potential
        A♭ = sv(A¹),
        A♯ = sv(Avec),
        DA♭ = sm(DA),
        DDA♭ = s3(DDA),
        φ = Num(φ⁰),

        # magnetic field
        B = Num(Babs),
        DB = sv(DBabs),
        DDB = sm(DDBabs),
        B♭ = sv(B¹),
        B♯ = sv(Bvec),
        B♮ = sv(Bphys),
        B♭♭ = sm(B²),
        DB♭ = sm(DB),

        # unit magnetic field
        b♭ = sv(b¹),
        b♯ = sv(bvec),
        b♮ = sv(bphys),
        Db♭ = sm(Db),
        Db♮ = sm(Dbphys),
        DDb♭ = s3(DDb),

        # perpendicular frame
        a♭ = sv(a¹),
        a♯ = sv(avec),
        a♮ = sv(aphys),
        c♭ = sv(c¹),
        c♯ = sv(cvec),
        c♮ = sv(cphys),

        # electric field
        E♭ = sv(E¹),
        E♯ = sv(Evec),
        DE♭ = sm(DE)
    )
end
