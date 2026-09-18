
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

        f(field, t, ξ)      # warm up
        @test @allocated(f(field, t, ξ)) == 0
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
# `Singular` diverges on the axis and its frame is left out, as it was before the rewrite.
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
    # identified by a hash of its body, so two built independently from the same expression are
    # `===` whether or not either came from the cache — which is exactly why the cache can
    # survive precompilation. It is still evidence that the parameters are arguments: baked in as
    # literals they would give `a` and `b` different bodies, and so different objects.
    #
    # Whether the cache was used is therefore observed through the cache itself.
    clear_field_cache!()
    @test isempty(ElectromagneticFields.FIELD_CACHE)

    uncached = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 1.7);
        cache = false)
    @test isempty(ElectromagneticFields.FIELD_CACHE)
    @test B♭(uncached, t, ξ) == B♭(b, t, ξ)

    # a cached build fills it, and a field rebuilt after a clear agrees with the one it replaces
    rebuilt = FieldFunctions(AxisymmetricTokamakCylindricalEquilibrium(6.2, 5.3, 1.7))
    @test !isempty(ElectromagneticFields.FIELD_CACHE)
    for name in ElectromagneticFields.FIELD_FUNCTION_NAMES
        f = getfield(ElectromagneticFields, name)
        @test f(rebuilt, t, ξ) == f(b, t, ξ)
    end
end
