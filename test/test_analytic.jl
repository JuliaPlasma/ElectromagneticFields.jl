
using ElectromagneticFields
using LinearAlgebra
using StaticArrays
using Test

using ElectromagneticFields: FIELD_FUNCTION_NAMES

# testing parameters
const t = 1.0
const ξ = [1.05, 0.5, 0.5]

generics() = (getfield(ElectromagneticFields, name) for name in FIELD_FUNCTION_NAMES)

"""
Allocations of one call to `f`, measured behind a function barrier.

`f` arrives untyped from `generics()`, so the call is dynamically dispatched and `@allocated`
applied to it directly charges for boxing the result rather than for the work — the barrier
`docs/src/interface.md` prescribes for users is needed here for the same reason.
"""
@noinline function allocations(f, field, t, ξ)
    f(field, t, ξ)      # warm up
    @allocated(f(field, t, ξ))
end

"""
Every generic accepts a coordinate vector and three scalars and gives the same answer; every one
returns the coordinates' float type; and none of them allocates.

The type check is what the conversion in `FieldFunction` exists for. A structurally constant body
such as `g♭` of a cartesian chart is emitted with `Int` literals, and left alone it would make the
type of a quantity depend on which equilibrium it came from.

The allocation check is what forces the coordinates through an `SVector`: `build_function` builds
its output container `similarto` its argument, so a plain `Vector` would return heap arrays and
allocate up to 6 KiB for a rank-three tensor.
"""
function test_interface(field, t, ξ)
    for f in generics()
        @test f(field, t, ξ...) == f(field, t, ξ)

        v = f(field, t, ξ)
        @test (v isa Number ? typeof(v) : eltype(v)) === eltype(ξ)

        @test allocations(f, field, t, ξ) == 0
    end

    # the domain bounds take no time
    @test rangemin(field, ξ) == rangemin(field, t, ξ)
    @test rangemax(field, ξ) == rangemax(field, t, ξ)
    @test rangemin(field, ξ...) == rangemin(field, t, ξ)
    @test rangemax(field, ξ...) == rangemax(field, t, ξ)
end

"The chart: the metric, the Jacobian, the volume element and the sign it throws away."
function test_chart(field, equ, t, ξ)
    @test from_cartesian(field, t, to_cartesian(field, t, ξ)) ≈ ξ atol = 1E-14

    G = g♭(field, t, ξ)
    Ḡ = g♯(field, t, ξ)
    F = DF(field, t, ξ)
    F̄ = DF̄(field, t, ξ)

    @test J(field, t, ξ) ≈ sqrt(det(F' * F)) atol = 1E-12

    # `J` is the unsigned volume element |det DF|, asserted just above. `orientation` carries the
    # sign that `J` throws away, and the two together must reproduce the signed determinant —
    # otherwise the Hodge star and the cross product, which are handed `orientation(equ) * J`, are
    # working in the wrong-handed frame.
    @test det(F) ≈ orientation(field) * J(field, t, ξ) atol = 1E-12
    @test orientation(field) ∈ (-1, +1)
    @test orientation(field) == ElectromagneticFields.orientation(equ)

    @test Ḡ ≈ inv(G) atol = 1E-12
    @test F̄ ≈ inv(F) atol = 1E-12
    @test F' * F ≈ G atol = 1E-12
    @test F * F̄ ≈ I atol = 1E-12
    @test F̄ * F̄' ≈ Ḡ atol = 1E-12
end

"Raising and lowering indices, and the three representations of the same field."
function test_representations(field, t, ξ)
    G = g♭(field, t, ξ)
    Ḡ = g♯(field, t, ξ)
    F = DF(field, t, ξ)
    F̄ = DF̄(field, t, ξ)

    for (cov, con) in ((A♭, A♯), (B♭, B♯), (E♭, E♯))
        @test con(field, t, ξ) ≈ Ḡ * cov(field, t, ξ) atol = 1E-12
        @test cov(field, t, ξ) ≈ G * con(field, t, ξ) atol = 1E-12
    end

    @test B♮(field, t, ξ) ≈ F̄' * B♭(field, t, ξ) atol = 1E-12
    @test B♮(field, t, ξ) ≈ F * B♯(field, t, ξ) atol = 1E-12

    @test B(field, t, ξ) ≈ sqrt(dot(B♯(field, t, ξ), B♭(field, t, ξ))) atol = 1E-12
    @test b♭(field, t, ξ) ≈ B♭(field, t, ξ) / B(field, t, ξ) atol = 1E-12
    @test b♯(field, t, ξ) ≈ B♯(field, t, ξ) / B(field, t, ξ) atol = 1E-12
    @test b♮(field, t, ξ) ≈ B♮(field, t, ξ) / B(field, t, ξ) atol = 1E-12

    # the two-form is antisymmetric by construction
    @test B♭♭(field, t, ξ) ≈ -transpose(B♭♭(field, t, ξ)) atol = 1E-14
end

"The perpendicular frame (a, b, c) is orthonormal in all three representations."
function test_frame(field, t, ξ)
    G = g♭(field, t, ξ)
    Ḡ = g♯(field, t, ξ)
    F = DF(field, t, ξ)
    F̄ = DF̄(field, t, ξ)

    for (cov, con, phys) in ((a♭, a♯, a♮), (b♭, b♯, b♮), (c♭, c♯, c♮))
        u, v, w = cov(field, t, ξ), con(field, t, ξ), phys(field, t, ξ)

        @test u ≈ G * v atol = 1E-14
        @test v ≈ Ḡ * u atol = 1E-14
        @test w ≈ F * v atol = 1E-14
        @test w ≈ F̄' * u atol = 1E-14

        @test dot(v, u) ≈ 1 atol = 1E-14
        @test dot(w, w) ≈ 1 atol = 1E-14
        @test u' * Ḡ * u ≈ 1 atol = 1E-14
        @test v' * G * v ≈ 1 atol = 1E-14
    end

    for (p, q) in ((a♮, b♮), (b♮, c♮), (c♮, a♮))
        @test dot(p(field, t, ξ), q(field, t, ξ)) ≈ 0 atol = 1E-14
    end
    for (p, q) in ((a♭, b♭), (b♭, c♭), (c♭, a♭))
        @test p(field, t, ξ)' * Ḡ * q(field, t, ξ) ≈ 0 atol = 1E-14
    end
    for (p, q) in ((a♯, b♯), (b♯, c♯), (c♯, a♯))
        @test p(field, t, ξ)' * G * q(field, t, ξ) ≈ 0 atol = 1E-14
    end
end

"""
`B` really is the curl of `A`: finite-difference ∇ × A in cartesian coordinates and compare with
the physical components of the generated field. This is the check that catches an orientation
error, which reverses `B` with no other visible symptom.
"""
function test_curl(field, t, ξ; h = 1E-5)
    x = to_cartesian(field, t, ξ)

    function A_cartesian(t, x)
        η = from_cartesian(field, t, x)
        F̄ = DF̄(field, t, η)
        F̄' * A♭(field, t, η)
    end

    ê(i) = SVector{3}(k == i ? one(eltype(x)) : zero(eltype(x)) for k in 1:3)
    ∂(i, j) = (A_cartesian(t, x .+ h .* ê(i))[j] - A_cartesian(t, x .- h .* ê(i))[j]) / (2h)

    curlA = [∂(2, 3) - ∂(3, 2), ∂(3, 1) - ∂(1, 3), ∂(1, 2) - ∂(2, 1)]
    Bcar = B♮(field, t, ξ)

    # central differences on an O(1) field carry an O(h²) truncation error; the tolerance is
    # scaled by |B| so that it means the same thing for the Dipole as for the ThetaPinch
    @test norm(curlA - Bcar) ≤ 1E-6 * max(norm(Bcar), 1)
end

# equilibrium, sample point, rangemin, rangemax, and whether the perpendicular frame is checked.
# `Singular` diverges on the axis, so its frame is left out.
const EQUILIBRIA = [
    ("ABC", ABCEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("AxisymmetricTokamakCartesian", AxisymmetricTokamakCartesianEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("AxisymmetricTokamakCylindrical", AxisymmetricTokamakCylindricalEquilibrium(), ξ,
        [-Inf, -Inf, 0.0], [+Inf, +Inf, 2π], true),
    ("AxisymmetricTokamakToroidal", AxisymmetricTokamakToroidalEquilibrium(), ξ,
        [-Inf, 0.0, 0.0], [+Inf, 2π, 2π], true),
    ("AxisymmetricTokamakToroidalRegularization",
        AxisymmetricTokamakToroidalRegularizationEquilibrium(), ξ,
        [-Inf, 0.0, 0.0], [+Inf, 2π, 2π], true),
    ("Dipole", DipoleField(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("PenningTrapUniform", PenningTrapUniformEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("PenningTrapBottle", PenningTrapBottleEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("PenningTrapAsymmetric", PenningTrapAsymmetricEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("QuadraticPotentials", QuadraticPotentialsField(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("Singular", SingularEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], false),
    ("SymmetricQuadratic", SymmetricQuadraticEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("ThetaPinch", ThetaPinchEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true),
    ("SolovevFRC", SolovevEquilibriumFRC(), ξ,
        [-Inf, -Inf, 0.0], [+Inf, +Inf, 2π], true),
    ("SolovevITER", SolovevEquilibriumITER(), ξ,
        [-Inf, -Inf, 0.0], [+Inf, +Inf, 2π], true),
    ("SolovevITERwXpoint", SolovevXpointEquilibriumITER(), ξ,
        [-Inf, -Inf, 0.0], [+Inf, +Inf, 2π], true),
    ("SolovevNSTX", SolovevEquilibriumNSTX(), ξ,
        [-Inf, -Inf, 0.0], [+Inf, +Inf, 2π], true),
    ("SolovevNSTXwXpoint", SolovevXpointEquilibriumNSTX(), ξ,
        [-Inf, -Inf, 0.0], [+Inf, +Inf, 2π], true),
    ("SolovevNSTXwDoubleXpoint", SolovevDoubleXpointEquilibriumNSTX(), ξ,
        [-Inf, -Inf, 0.0], [+Inf, +Inf, 2π], true),
    ("SolovevSymmetric", SolovevSymmetricEquilibrium(), ξ,
        [-Inf, -Inf, -Inf], [+Inf, +Inf, +Inf], true)
]

const FIELDS = Dict{String, Any}()

for (name, equ, p, rmin, rmax, frame) in EQUILIBRIA
    field = FieldFunctions(equ)
    FIELDS[name] = field

    @testset "$(rpad(name, 60))" begin
        test_interface(field, t, p)
        test_chart(field, equ, t, p)
        test_representations(field, t, p)
        frame && test_frame(field, t, p)
        test_curl(field, t, p)

        @test rangemin(field, t, p) == rmin
        @test rangemax(field, t, p) == rmax
    end
end

# test correctness of some of the magnetic fields

function test_axisymmetric_tokamak_cartesian_equilibrium(
        field, t = 0.0, x = [1.5, 0.0, 0.5])
    par = parameters(field)
    crd = coordinates(field)

    @test B♯(field, t, x) ≈ B♭(field, t, x) atol = 1E-16

    @test B♭(field, t, x)[1] ≈
          -par.B₀ / par.q₀ * (par.q₀ * par.R₀ * crd.Y(t, x) + crd.X(t, x) * crd.Z(t, x)) /
          crd.R(t, x)^2 atol = 1E-16
    @test B♭(field, t, x)[2] ≈
          +par.B₀ / par.q₀ * (par.q₀ * par.R₀ * crd.X(t, x) - crd.Y(t, x) * crd.Z(t, x)) /
          crd.R(t, x)^2 atol = 1E-16
    @test B♭(field, t, x)[3] ≈
          +par.B₀ / par.q₀ * (crd.R(t, x) - par.R₀) / crd.R(t, x) atol = 1E-16
end

function test_axisymmetric_tokamak_cylindrical_equilibrium(
        field, t = 0.0, x = [1.5, 0.5, π / 5])
    par = parameters(field)
    crd = coordinates(field)

    @test B♯(field, t, x)[1] == -par.B₀ / par.q₀ * crd.Z(t, x) / crd.R(t, x)
    @test B♯(field, t, x)[2] == +par.B₀ / par.q₀ * (crd.R(t, x) - par.R₀) / crd.R(t, x)
    @test B♯(field, t, x)[3] == +par.B₀ * par.R₀ / crd.R(t, x)^2

    @test B♭(field, t, x)[1] == -par.B₀ / par.q₀ * crd.Z(t, x) / crd.R(t, x)
    @test B♭(field, t, x)[2] == +par.B₀ / par.q₀ * (crd.R(t, x) - par.R₀) / crd.R(t, x)
    @test B♭(field, t, x)[3] == +par.B₀ * par.R₀
end

function test_axisymmetric_tokamak_toroidal_equilibrium(
        field, t = 0.0, x = [0.5, π / 10, π / 5])
    par = parameters(field)
    crd = coordinates(field)

    @test B♯(field, t, x)[1] == 0
    @test B♯(field, t, x)[2] == +par.B₀ / par.q₀ / crd.R(t, x)
    @test B♯(field, t, x)[3] ≈ +par.B₀ * par.R₀ / crd.R(t, x)^2 atol = 1E-14

    @test B♭(field, t, x)[1] == 0
    @test B♭(field, t, x)[2] == +par.B₀ / par.q₀ * crd.r(t, x)^2 / crd.R(t, x)
    @test B♭(field, t, x)[3] ≈ +par.B₀ * par.R₀ atol = 1E-14
end

"""
The regularised chart carries the same magnetic field as
`AxisymmetricTokamakToroidalEquilibrium`, in a gauge whose poloidal vector potential is regular on
the magnetic axis, so it must reproduce that chart's field values exactly. Only the tolerances
differ: the `1/cos²θ` gauge makes the generated expressions less well conditioned, so these are
approximate where the unregularised ones are exact.
"""
function test_axisymmetric_tokamak_toroidal_regularization_equilibrium(
        field, t = 0.0, x = [0.5, π / 10, π / 5])
    par = parameters(field)
    crd = coordinates(field)

    @test B♯(field, t, x)[1] ≈ 0 atol = 1E-14
    @test B♯(field, t, x)[2] ≈ +par.B₀ / par.q₀ / crd.R(t, x) atol = 1E-14
    @test B♯(field, t, x)[3] ≈ +par.B₀ * par.R₀ / crd.R(t, x)^2 atol = 1E-14

    @test B♭(field, t, x)[1] ≈ 0 atol = 1E-14
    @test B♭(field, t, x)[2] ≈ +par.B₀ / par.q₀ * crd.r(t, x)^2 / crd.R(t, x) atol = 1E-14
    @test B♭(field, t, x)[3] ≈ +par.B₀ * par.R₀ atol = 1E-14
end

function test_symmetric_quadratic_equilibrium(field, t = 0.0, x = [1.0, 0.5, 0.5])
    par = parameters(field)
    crd = coordinates(field)

    @test B♯(field, t, x) == B♭(field, t, x)

    @test B♭(field, t, x)[1] == 0
    @test B♭(field, t, x)[2] == 0
    @test B♭(field, t, x)[3] == par.B₀ * (1 + crd.X(t, x)^2 + crd.Y(t, x)^2)

    @test B(field, t, x) == par.B₀ * (1 + crd.X(t, x)^2 + crd.Y(t, x)^2)

    @test b♯(field, t, x) == [0, 0, 1]
    @test b♭(field, t, x) == [0, 0, 1]
end

function test_theta_pinch_equilibrium(field, t = 0.0, x = [1.0, 0.5, 0.5])
    par = parameters(field)

    @test B♯(field, t, x) == B♭(field, t, x)
    @test B♭(field, t, x) == [0, 0, par.B₀]
    @test B(field, t, x) == par.B₀
    @test b♯(field, t, x) == [0, 0, 1]
    @test b♭(field, t, x) == [0, 0, 1]
end

function test_abc_equilibrium(field, t = 0.0, x = [1.0, 0.5, 0.5])
    par = parameters(field)
    crd = coordinates(field)

    @test B♯(field, t, x) == B♭(field, t, x)
    @test B♭(field, t, x) == A♭(field, t, x)

    @test B♭(field, t, x)[1] == par.a₀ * sin(crd.Z(t, x)) + par.c₀ * cos(crd.Y(t, x))
    @test B♭(field, t, x)[2] == par.b₀ * sin(crd.X(t, x)) + par.a₀ * cos(crd.Z(t, x))
    @test B♭(field, t, x)[3] == par.c₀ * sin(crd.Y(t, x)) + par.b₀ * cos(crd.X(t, x))
end

@testset "$(rpad("Magnetic Fields", 60))" begin
    test_axisymmetric_tokamak_cartesian_equilibrium(FIELDS["AxisymmetricTokamakCartesian"])
    test_axisymmetric_tokamak_cylindrical_equilibrium(FIELDS["AxisymmetricTokamakCylindrical"])
    test_axisymmetric_tokamak_toroidal_equilibrium(FIELDS["AxisymmetricTokamakToroidal"])
    test_axisymmetric_tokamak_toroidal_regularization_equilibrium(FIELDS["AxisymmetricTokamakToroidalRegularization"])
    test_symmetric_quadratic_equilibrium(FIELDS["SymmetricQuadratic"])
    test_theta_pinch_equilibrium(FIELDS["ThetaPinch"])
    test_abc_equilibrium(FIELDS["ABC"])
end

"""
Two charts of the same physical field, compared at the same physical point: the contravariant
components transform with `DF`, the covariant ones with `DF̄'`, and their contraction is a scalar
and so chart-independent.
"""
function test_chart_consistency(field, field_car, t, p)
    x = to_cartesian(field, t, p)

    @test dot(B♯(field, t, p), B♭(field, t, p)) ≈
          dot(B♯(field_car, t, x), B♭(field_car, t, x)) atol = 1E-12
    @test DF(field, t, p) * B♯(field, t, p) ≈ B♯(field_car, t, x) atol = 1E-12
    @test DF̄(field, t, p)' * B♭(field, t, p) ≈ B♭(field_car, t, x) atol = 1E-12
end

@testset "$(rpad("Consistency", 60))" begin
    car = FIELDS["AxisymmetricTokamakCartesian"]
    test_chart_consistency(FIELDS["AxisymmetricTokamakCylindrical"], car, 0.0,
        [1.5, 0.5, π / 5])
    test_chart_consistency(FIELDS["AxisymmetricTokamakToroidal"], car, 0.0,
        [0.5, π / 10, π / 5])
    # the regularised chart shares the toroidal chart's coordinates, so the same check applies
    test_chart_consistency(FIELDS["AxisymmetricTokamakToroidalRegularization"], car, 0.0,
        [0.5, π / 10, π / 5])
end

# `A♭[3]` is the poloidal flux function of the axisymmetric equilibria, and it is what the plotting
# extension contours. The defining property is that the magnetic field lies in its level surfaces,
# `B · ∇A₃ = 0`. The second half of each case guards the distinction that makes this worth
# asserting: the physical toroidal component `A₃ / R` is a plausible-looking stand-in that does not
# have the property. `R` there is the major radius in each chart's own coordinates, not `ξ₁` — in
# the toroidal chart `ξ₁` is `r`, and dividing by it leaves a flux label behind.

@testset "$(rpad("A₃ is a flux label for the axisymmetric equilibria", 60))" begin
    for (field, p) in (
        (FIELDS["AxisymmetricTokamakCylindrical"], [1.1, 0.2, 0.3]),
        (FIELDS["AxisymmetricTokamakToroidal"], [0.2, 0.7, 0.3]),
        (FIELDS["SolovevITER"], [1.1, 0.2, 0.3])
    )
        Bcon = B♯(field, t, p)

        # central differences, so the tolerance is set by the truncation error rather than by ε
        h = 1E-6
        ê(i) = [k == i ? h : zero(h) for k in 1:3]
        ∇(f) = [(f(p .+ ê(i)) - f(p .- ê(i))) / 2h for i in 1:3]

        ψ(q) = A♭(field, t, q)[3]
        @test dot(Bcon, ∇(ψ)) ≈ 0 atol = 1E-8

        # the quantity that is *not* a flux label, at a point where the difference shows
        physical(q) = A♭(field, t, q)[3] / coordinates(field).R(t, q)
        @test !isapprox(dot(Bcon, ∇(physical)), 0; atol = 1E-8)
    end
end

# The Solov'ev constructors solve a dense linear system for the coefficients `c`, and part of what
# that system says is that ψ vanishes at the boundary points of the Cerfon-Freidberg cross-section
# model. Those are the rows checked below, and nothing else in this file constrains `c` at all:
# every identity tested above holds for any `c` whatever.
#
# The cover is partial, and cannot be otherwise. The system is square and solved exactly, so `c`
# satisfies whichever rows it is given; a fault in one of the rows that constrain a first or second
# derivative moves `c` without moving these residuals. That is three of seven rows here for the
# symmetric and double-X families, and four of twelve for the X-point family.
#
# The three families impose different points, which is why they are listed per case rather than
# taken from the type — `SolovevDoubleXpointEquilibrium` returns a `SolovevXpointEquilibrium` too,
# and puts the X-point where the other two put the upper boundary point.
#
# `1E-13` is set by the measurement rather than by taste: the residuals run to 5E-16, against a ψ
# that measures 2E-3 at a sample point for the two ITER cases and 2E-1 to 4E-1 for the other four,
# and evaluating one of these at a point the system does not constrain gives 5E-3.

@testset "$(rpad("The Solov'ev boundary conditions are satisfied", 60))" begin
    outer(p) = [1 + p.ϵ, 0.0, 0.0]
    inner(p) = [1 - p.ϵ, 0.0, 0.0]
    upper(p) = [1 - p.δ * p.ϵ, p.κ * p.ϵ, 0.0]
    xpoint(p) = [p.xsep, p.ysep, 0.0]

    for (name, points) in (
        ("SolovevITER", (outer, inner, upper)),
        ("SolovevNSTX", (outer, inner, upper)),
        ("SolovevFRC", (outer, inner, upper)),
        ("SolovevITERwXpoint", (outer, inner, upper, xpoint)),
        ("SolovevNSTXwXpoint", (outer, inner, upper, xpoint)),
        ("SolovevNSTXwDoubleXpoint", (outer, inner, xpoint))
    )
        field = FIELDS[name]
        p = parameters(field)
        for point in points
            @test A♭(field, t, point(p))[3] ≈ 0 atol = 1E-13
        end
    end
end

# The frame is built from the basis vector ∂₁ unless the trace shows ∂₁ × b to be zero, and the
# choice is not revisited at each point. So it reverses on a path through the points where b ∥ ∂₁,
# and no choice of basis vector avoids that for a field whose b takes every direction. The dipole
# is one: b = ∓∂ₓ on the two lines y = 0, x = ±√2 z. In the cartesian chart
# a = ∂ₓ × b / |∂ₓ × b|, which at y = 0 is (0, -sign(B_z), 0). On x = +√2 z, where b ≈ -∂ₓ,
# c = b × a is then (0, 0, sign(B_z)) up to the small b_z. The expected sign comes from the closed
# form of B_z in the `DipoleField` docstring rather than from the generated code.

@testset "$(rpad("The dipole frame turns across b ∥ ∂₁", 60))" begin
    field = FIELDS["Dipole"]
    B₀ = parameters(field).B₀
    Bz(q) = -B₀ * (2q[3]^2 - q[1]^2 - q[2]^2) / norm(q)^5

    below = [sqrt(2) - 1E-8, 0.0, 1.0]
    above = [sqrt(2) + 1E-8, 0.0, 1.0]

    for q in (below, above)
        s = sign(Bz(q))
        @test a♯(field, t, q) ≈ [0, -s, 0]
        @test c♯(field, t, q) ≈ [0, 0, s] atol = 1E-6
    end

    @test a♯(field, t, below) ≈ -a♯(field, t, above)
    @test c♯(field, t, below) ≈ -c♯(field, t, above) atol = 1E-6
end

# `periodic` is answered per chart family rather than per equilibrium, because it is a property of
# the chart, and it has no fallback. So there are three things to check: that each family gives the
# right answer, that the families between them account for every shipped equilibrium, and that a
# chart which has not answered fails loudly instead of reporting no periodicity.
#
# A chart that has not answered is the case worth a test of its own. An all-`false` default would
# be indistinguishable, at every call site, from a chart that really has no periodic coordinate.

module NoChartAnswer

using ElectromagneticFields: AnalyticEquilibrium
import ElectromagneticFields: A₁, A₂, A₃, x¹, x², x³, ξ¹, ξ², ξ³, J

# A chart of its own — subtyping `AnalyticEquilibrium` rather than `CartesianEquilibrium` is what
# keeps it away from the cartesian family's method. It is the identity map, so it is a well-formed
# chart in every respect but the one being tested.
struct UnansweredChartEquilibrium{T <: Number} <: AnalyticEquilibrium
    name::String
    B₀::T

    function UnansweredChartEquilibrium{T}(B₀::T) where {T <: Number}
        new("UnansweredChartEquilibrium", B₀)
    end
end

function UnansweredChartEquilibrium(B₀::T) where {T <: Number}
    UnansweredChartEquilibrium{T}(B₀)
end

x¹(ξ::AbstractVector, ::UnansweredChartEquilibrium) = ξ[1]
x²(ξ::AbstractVector, ::UnansweredChartEquilibrium) = ξ[2]
x³(ξ::AbstractVector, ::UnansweredChartEquilibrium) = ξ[3]
ξ¹(x::AbstractVector, ::UnansweredChartEquilibrium) = x[1]
ξ²(x::AbstractVector, ::UnansweredChartEquilibrium) = x[2]
ξ³(x::AbstractVector, ::UnansweredChartEquilibrium) = x[3]
J(x::AbstractVector, ::UnansweredChartEquilibrium) = one(eltype(x))

A₁(x::AbstractVector, equ::UnansweredChartEquilibrium) = -equ.B₀ * x[2] / 2
A₂(x::AbstractVector, equ::UnansweredChartEquilibrium) = +equ.B₀ * x[1] / 2
A₃(x::AbstractVector, ::UnansweredChartEquilibrium) = zero(eltype(x))

end

@testset "$(rpad("Periodicity is a property of the chart", 60))" begin
    cartesian = ["ABC", "AxisymmetricTokamakCartesian", "Dipole", "PenningTrapUniform",
        "PenningTrapBottle", "PenningTrapAsymmetric", "QuadraticPotentials", "Singular",
        "SymmetricQuadratic", "ThetaPinch", "SolovevSymmetric"]
    cylindrical = ["AxisymmetricTokamakCylindrical", "SolovevFRC", "SolovevITER",
        "SolovevITERwXpoint", "SolovevNSTX", "SolovevNSTXwXpoint", "SolovevNSTXwDoubleXpoint"]
    toroidal = ["AxisymmetricTokamakToroidal", "AxisymmetricTokamakToroidalRegularization"]

    # the three families account for every equilibrium the table above builds, so a new one cannot
    # be added without landing in a family here or failing this line
    @test sort(vcat(cartesian, cylindrical, toroidal)) == sort([e[1] for e in EQUILIBRIA])

    for (names, expected) in ((cartesian, SVector(false, false, false)),
        (cylindrical, SVector(false, false, true)),
        (toroidal, SVector(false, true, true)))
        for name in names
            @test periodic(FIELDS[name]) == expected
        end
    end

    # the field stores what the equilibrium says, with the element type it says it in
    for (name, equ, _, _, _, _) in EQUILIBRIA
        @test periodic(FIELDS[name]) === periodic(equ)
        @test periodic(equ) isa SVector{3, Bool}
    end

    # no fallback: a chart nobody has answered for raises rather than answering `false` everywhere
    equ = NoChartAnswer.UnansweredChartEquilibrium(2.0)
    @test_throws MethodError periodic(equ)

    # and building a field from it raises too. `periodic` is the last thing the constructor
    # calls, after the whole symbolic trace, so the error is matched on its function as well as on
    # its type: a bare `MethodError` would equally match one raised earlier by some other part of
    # the chart interface, and the test would pass without reaching the line it is about.
    err = try
        FieldFunctions(equ)
        nothing
    catch e
        e
    end
    @test err isa MethodError && err.f === periodic
end

# A perturbation is combined with its equilibrium symbolically, before any code is generated, so
# the perturbed field is a `FieldFunctions` like any other. `EzCosZPerturbation` contributes only a
# scalar potential, so `B` is untouched and `E` is not.

@testset "$(rpad("Perturbation", 60))" begin
    equ = ThetaPinchEquilibrium()
    pert = EzCosZPerturbation(2.0)

    plain = FIELDS["ThetaPinch"]
    perturbed = FieldFunctions(equ, pert)

    @test perturbation(perturbed) === pert
    @test equilibrium(perturbed) === equ

    @test B♭(perturbed, t, ξ) ≈ B♭(plain, t, ξ)
    @test φ(plain, t, ξ) == 0
    @test φ(perturbed, t, ξ) ≈ 2.0 / (2π) * sin(2π * ξ[3])
    @test E♭(perturbed, t, ξ)[3] ≈ -2.0 * cos(2π * ξ[3])
    @test E♭(plain, t, ξ) == [0, 0, 0]
end

# The generated code takes the equilibrium's parameters as an argument instead of having them
# baked in, which is what lets one compiled function serve every parameter value of a type — and
# in turn what makes the cache, and the precompilation that fills it, correct.

# How far two builds of one equilibrium may drift apart. Reassociating the same arithmetic moves
# the last bit or two, and how far depends on what the platform contracts; a wrong formula moves it
# by O(1), so this is tight enough to catch one and four orders tighter than `≈` on its own.
const REBUILD_RTOL = 1.0e-12

@testset "$(rpad("Symbolic parameters and the cache", 60))" begin
    a = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(1.0, 1.0, 2.0))
    b = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 1.7))
    c = FieldFunctions(AxisymmetricTokamakToroidalEquilibrium())

    # same type, different parameters: the very same generated function, different answers
    for name in ElectromagneticFields.FIELD_FUNCTION_NAMES
        @test functions(a)[name].f === functions(b)[name].f
    end
    @test B♭(a, t, ξ) != B♭(b, t, ξ)

    # a different equilibrium type must not share it
    @test functions(c).B♭.f !== functions(a).B♭.f

    # the parameters are the equilibrium's own, and the values travel with the field
    @test parameters(b) == (R₀ = 6.2, B₀ = 5.3, q₀ = 1.7)
    @test functions(b).B♭.p == SVector(6.2, 5.3, 1.7)

    # `c` is derived rather than chosen, but `A₃` reads it, so it is a parameter like the rest
    @test :c ∈ keys(parameters(FIELDS["SolovevITER"]))

    # Note what the identity above does and does not witness. A `RuntimeGeneratedFunction` is
    # identified by a hash of its body, so two built from the same expression are `===` whether or
    # not either came from the cache — which is exactly why the cache can survive precompilation.
    # It is still evidence that the parameters are arguments: baked in as literals they would give
    # `a` and `b` different bodies, and so different objects.
    #
    # It says nothing about two separate traces of the same equilibrium, which need not produce the
    # same expression at all: the simplifier is free to choose any equivalent form, and
    # SymbolicUtils 4.46.8 chooses a different one for about half of these functions. So a rebuild
    # guarantees the same field, not the same rounding — see the comparison below.
    #
    # Whether the cache was used is therefore observed through the cache itself.
    clear_field_cache!()
    @test isempty(ElectromagneticFields.FIELD_CACHE)

    uncached = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 1.7);
        cache = false)
    @test isempty(ElectromagneticFields.FIELD_CACHE)
    @test B♭(uncached, t, ξ) ≈ B♭(b, t, ξ) rtol=REBUILD_RTOL

    # a cached build fills it, and a field rebuilt after a clear agrees with the one it replaces.
    # The parameters are compared exactly, because those are the equilibrium's own numbers and a
    # rebuild copies them. The values are compared to a few ULP, because the two bodies may be
    # different orderings of the same arithmetic; a wrong formula is off by O(1) and still fails.
    rebuilt = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 1.7))
    @test !isempty(ElectromagneticFields.FIELD_CACHE)
    for name in ElectromagneticFields.FIELD_FUNCTION_NAMES
        @test functions(rebuilt)[name].p == functions(b)[name].p
        f = getfield(ElectromagneticFields, name)
        @test f(rebuilt, t, ξ) ≈ f(b, t, ξ) rtol=REBUILD_RTOL
    end
end

# A type with two vector parameters is the case a cache key over the flattened total cannot
# separate. Lengths `(2, 3)` and `(3, 2)` both flatten to five slots, and the generated code reads
# those slots positionally, so the second field would be served the first one's code and return
# numbers computed from the wrong parameters — no bounds error, no `MethodError`, and `parameters`
# still showing the right struct. Nothing shipped here can reach it: the two structs that carry a
# vector parameter carry exactly one each, so their total fixes its length and with it the split.
module TwoVectorField

using ElectromagneticFields: CartesianEquilibrium, CartesianPerturbation, X, Y, Z
import ElectromagneticFields: A₁, A₂, A₃, φ

struct TwoVectorEquilibrium{T <: Number} <: CartesianEquilibrium
    name::String
    a::Vector{T}
    b::Vector{T}

    function TwoVectorEquilibrium{T}(a::Vector{T}, b::Vector{T}) where {T <: Number}
        new("TwoVectorEquilibrium", a, b)
    end
end

function TwoVectorEquilibrium(a::Vector{T}, b::Vector{T}) where {T <: Number}
    TwoVectorEquilibrium{T}(a, b)
end

A₁(x::AbstractVector, ::TwoVectorEquilibrium) = zero(eltype(x))
A₂(x::AbstractVector, ::TwoVectorEquilibrium) = zero(eltype(x))

# reads both parameters and weights them differently, so reading the slots under the other split
# changes the answer rather than merely reordering it
function A₃(x::AbstractVector, equ::TwoVectorEquilibrium)
    sum(equ.a) * X(x, equ) + 2 * sum(equ.b) * Y(x, equ)
end

# One vector parameter each, so the pair is what varies the split between the equilibrium's
# parameters and the perturbation's. `pvalues` concatenates the two, so a key over the total
# cannot see where the boundary falls.
struct OneVectorEquilibrium{T <: Number} <: CartesianEquilibrium
    name::String
    a::Vector{T}

    function OneVectorEquilibrium{T}(a::Vector{T}) where {T <: Number}
        new("OneVectorEquilibrium", a)
    end
end

OneVectorEquilibrium(a::Vector{T}) where {T <: Number} = OneVectorEquilibrium{T}(a)

A₁(x::AbstractVector, ::OneVectorEquilibrium) = zero(eltype(x))
A₂(x::AbstractVector, ::OneVectorEquilibrium) = zero(eltype(x))
A₃(x::AbstractVector, equ::OneVectorEquilibrium) = sum(equ.a) * X(x, equ)

struct OneVectorPerturbation{T <: Number} <: CartesianPerturbation
    name::String
    e::Vector{T}

    function OneVectorPerturbation{T}(e::Vector{T}) where {T <: Number}
        new("OneVectorPerturbation", e)
    end
end

OneVectorPerturbation(e::Vector{T}) where {T <: Number} = OneVectorPerturbation{T}(e)

φ(x::AbstractVector, pert::OneVectorPerturbation) = sum(pert.e) * Z(x, pert)

# a parameter field wide enough to hold either a scalar or a vector, which is what tells the
# encoding of a scalar apart from that of an empty vector. No field is built from it.
struct LooseEquilibrium <: CartesianEquilibrium
    name::String
    p::Any
end

end

@testset "$(rpad("The cache key carries the parameter shape", 60))" begin
    equ23 = TwoVectorField.TwoVectorEquilibrium([1.0, 2.0], [3.0, 4.0, 5.0])
    equ32 = TwoVectorField.TwoVectorEquilibrium([1.0, 2.0, 3.0], [4.0, 5.0])

    @test ElectromagneticFields.parameter_shape(equ23) == (2, 3)
    @test ElectromagneticFields.parameter_shape(equ32) == (3, 2)
    @test ElectromagneticFields.parameter_shape(ThetaPinchEquilibrium()) == (-1,)
    @test ElectromagneticFields.parameter_shape(ZeroPerturbation()) == ()

    # a scalar is not a vector of length zero. Were they one encoding, the shape would be coarser
    # than the total it replaces, and this pair would share a cache entry that the total separates.
    scalar = TwoVectorField.LooseEquilibrium("Loose", 1.0)
    empty = TwoVectorField.LooseEquilibrium("Loose", Float64[])
    @test ElectromagneticFields.parameter_shape(scalar) !=
          ElectromagneticFields.parameter_shape(empty)
    @test length(ElectromagneticFields.parameter_values(scalar)) !=
          length(ElectromagneticFields.parameter_values(empty))

    # the two shapes flatten to the same number of slots: a total alone cannot separate them
    @test length(ElectromagneticFields.parameter_values(equ23)) ==
          length(ElectromagneticFields.parameter_values(equ32))

    entries = length(ElectromagneticFields.FIELD_CACHE)
    f23 = FieldFunctions(equ23)
    f32 = FieldFunctions(equ32)

    # one entry each, with different code, rather than one entry serving both
    @test length(ElectromagneticFields.FIELD_CACHE) == entries + 2
    @test functions(f23).A♭.f !== functions(f32).A♭.f

    # so each field reads its own parameters. `A₃ = sum(a) x + 2 sum(b) y`, and the numbers are
    # small integers, so the comparison is exact however the trace associates the sums.
    ζ = [1.0, 1.0, 0.0]
    @test A♭(f23, t, ζ)[3] == sum(equ23.a) + 2 * sum(equ23.b)
    @test A♭(f32, t, ζ)[3] == sum(equ32.a) + 2 * sum(equ32.b)

    # the perturbation's shape is a key element of its own, so the boundary between the two
    # structs' parameters is visible. `pvalues` concatenates them, so both pairs below flatten to
    # the same five values and only the boundary moves.
    eq2 = TwoVectorField.OneVectorEquilibrium([1.0, 2.0])
    eq3 = TwoVectorField.OneVectorEquilibrium([1.0, 2.0, 3.0])
    pert3 = TwoVectorField.OneVectorPerturbation([3.0, 4.0, 5.0])
    pert2 = TwoVectorField.OneVectorPerturbation([4.0, 5.0])

    @test length(ElectromagneticFields.parameter_values(eq2)) +
          length(ElectromagneticFields.parameter_values(pert3)) ==
          length(ElectromagneticFields.parameter_values(eq3)) +
          length(ElectromagneticFields.parameter_values(pert2))

    split = length(ElectromagneticFields.FIELD_CACHE)
    g23 = FieldFunctions(eq2, pert3)
    g32 = FieldFunctions(eq3, pert2)
    @test length(ElectromagneticFields.FIELD_CACHE) == split + 2

    ζ₃ = [1.0, 1.0, 1.0]
    @test A♭(g23, t, ζ₃)[3] == sum(eq2.a)
    @test A♭(g32, t, ζ₃)[3] == sum(eq3.a)
    @test φ(g23, t, ζ₃) == sum(pert3.e)
    @test φ(g32, t, ζ₃) == sum(pert2.e)

    # the key is no finer than it must be: a second field of one shape is still a lookup
    filled = length(ElectromagneticFields.FIELD_CACHE)
    other23 = TwoVectorField.TwoVectorEquilibrium([6.0, 7.0], [8.0, 9.0, 10.0])
    again = FieldFunctions(other23)
    @test length(ElectromagneticFields.FIELD_CACHE) == filled
    @test functions(again).A♭.f === functions(f23).A♭.f
    @test A♭(again, t, ζ)[3] == sum(other23.a) + 2 * sum(other23.b)
end

# `get_parameters` is called on the instance, so one type can name a different set of parameters
# for each instance. The shape cannot see that: `(:a, :b)` and `(:b, :c)` are both `(-1, -1)`. The
# test type chooses its parameters per instance and leaves the third field frozen as a literal,
# which is what such a type does, and it needs its own `symbolic_copy`, because the default one
# passes the parameters to the constructor in order.

module ChosenParameterField

using ElectromagneticFields: CartesianEquilibrium, Symbolics, X, Y, Z
import ElectromagneticFields: A₁, A₂, A₃, get_parameters, symbolic_copy

struct ChosenEquilibrium{T <: Number} <: CartesianEquilibrium
    name::String
    chosen::NTuple{2, Symbol}
    a::T
    b::T
    c::T
end

get_parameters(equ::ChosenEquilibrium) = equ.chosen

function symbolic_copy(equ::ChosenEquilibrium, prefix::Symbol)
    symbols = map(name -> Symbolics.variable(Symbol(prefix, :_, name)), equ.chosen)
    member(name) = name in equ.chosen ? symbols[findfirst(==(name), equ.chosen)] :
                   Symbolics.Num(getfield(equ, name))
    copy = ChosenEquilibrium{Symbolics.Num}(
        equ.name, equ.chosen, member(:a), member(:b), member(:c))
    copy, collect(symbols)
end

A₁(x::AbstractVector, ::ChosenEquilibrium) = zero(eltype(x))
A₂(x::AbstractVector, ::ChosenEquilibrium) = zero(eltype(x))

# the three weights differ, so reading one parameter's slot as another's changes the answer
function A₃(x::AbstractVector, equ::ChosenEquilibrium)
    equ.a * X(x, equ) + 2 * equ.b * Y(x, equ) + 4 * equ.c * Z(x, equ)
end

end

@testset "$(rpad("The cache key carries the parameter names", 60))" begin
    ab = ChosenParameterField.ChosenEquilibrium("Chosen", (:a, :b), 1.0, 2.0, 3.0)
    bc = ChosenParameterField.ChosenEquilibrium("Chosen", (:b, :c), 1.0, 2.0, 3.0)

    # the same type and the same shape, so only the names tell the two apart
    @test typeof(ab) == typeof(bc)
    @test ElectromagneticFields.parameter_shape(ab) ==
          ElectromagneticFields.parameter_shape(bc)
    @test ElectromagneticFields.parameter_names(ab) !=
          ElectromagneticFields.parameter_names(bc)

    entries = length(ElectromagneticFields.FIELD_CACHE)
    fab = FieldFunctions(ab)
    fbc = FieldFunctions(bc)

    # one entry each, rather than `bc` served the code traced for `ab`
    @test length(ElectromagneticFields.FIELD_CACHE) == entries + 2
    @test functions(fab).A♭.f !== functions(fbc).A♭.f

    # `A₃ = a x + 2 b y + 4 c z`, with small integers, so the comparison is exact. Served `ab`'s
    # code, `bc` would read `b` as `a` and `c` as `b`, and give 2 + 6 + 12 = 20 here.
    ζ = [1.0, 1.0, 1.0]
    @test A♭(fab, t, ζ)[3] == 1 + 4 + 12
    @test A♭(fbc, t, ζ)[3] == 1 + 4 + 12

    # the key is no finer than it must be: the same names are still a lookup, whatever the values
    filled = length(ElectromagneticFields.FIELD_CACHE)
    other = ChosenParameterField.ChosenEquilibrium("Chosen", (:a, :b), 5.0, 6.0, 3.0)
    again = FieldFunctions(other)
    @test length(ElectromagneticFields.FIELD_CACHE) == filled
    @test functions(again).A♭.f === functions(fab).A♭.f
    @test A♭(again, t, ζ)[3] == 5 + 12 + 12
end
