
using Combinatorics
using LinearAlgebra
using StaticArrays
using Symbolics

abstract type AnalyticField <: ElectromagneticField end
abstract type AnalyticEquilibrium <: AnalyticField end
abstract type AnalyticPerturbation <: AnalyticField end

function get_functions end
function get_parameters end

# Coordinate helpers an equilibrium may define for its own chart, each as a method taking
# `(x, equ)`. They were module-local in every field module before the modules were flattened, so
# the same name now carries one method per equilibrium. `get_functions` lists which of them to
# generate into `coordinates(field)`. They are deliberately not exported — `r`, `θ` and `ϕ` would
# collide with practically any caller.
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

function GeometricBase.periodicity(x::AbstractVector{T}, ::AnalyticField) where {T}
    (-Inf * ones(T, 4), +Inf * ones(T, 4))
end

from_cartesian(x::AbstractVector, equ::AnalyticField) = [ξ¹(x, equ), ξ²(x, equ), ξ³(x, equ)]
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

"Returns the i-th component of the physical coordinate representation of the one-form α"
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
    generate_field_expressions(equ, pert)

Build the symbolic expression for every quantity [`FieldFunctions`](@ref) stores, returned as a
`NamedTuple` whose keys are that struct's field names.

Each value is a `Num`, or a `StaticArray` of `Num` of the quantity's tensor rank: `SVector{3}` for
a one-form or a vector, `SMatrix{3,3}` for a Jacobian or the metric, `SArray{Tuple{3,3,3}}` and
`SArray{Tuple{3,3,3,3}}` for the higher derivatives. The container matters — `build_function`
builds its output `similarto` what it is given.
"""
function generate_field_expressions(equ::AnalyticEquilibrium, pert::AnalyticPerturbation)
    # Symbols for time t and coordinates x = (x₁, x₂, x₃), ξ = (ξ₁, ξ₂, ξ₃).
    Symbolics.@variables t x₁ x₂ x₃ ξ₁ ξ₂ ξ₃
    x = [x₁, x₂, x₃]
    ξ = [ξ₁, ξ₂, ξ₃]

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
        if !all(iszero, Symbolics.simplify.(avec))
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
        arguments = (t, ξ),
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
