
# some convenience functions
function structname(equ)
    if occursin('.', equ)
        return equ[findlast(isequal('.'), equ)+1:end]
    else
        return equ
    end
end

teststring(equ) = structname(string(equ))
teststring(equ, pert) = structname(string(equ)) * " + " * structname(string(pert))


function test_equilibrium(equ_mod, t, ξ)
    @test equ_mod.X(t, ξ...) == equ_mod.X(t,ξ)
    @test equ_mod.Y(t, ξ...) == equ_mod.Y(t,ξ)
    @test equ_mod.Z(t, ξ...) == equ_mod.Z(t,ξ)

    @test equ_mod.J(t, ξ...) == equ_mod.J(t,ξ)
    @test equ_mod.B(t, ξ...) == equ_mod.B(t,ξ)
    @test equ_mod.φ(t, ξ...) == equ_mod.φ(t,ξ)

    @test equ_mod.A₁(t, ξ...) == equ_mod.A₁(t,ξ)
    @test equ_mod.A₂(t, ξ...) == equ_mod.A₂(t,ξ)
    @test equ_mod.A₃(t, ξ...) == equ_mod.A₃(t,ξ)

    @test equ_mod.B₁(t, ξ...) == equ_mod.B₁(t,ξ)
    @test equ_mod.B₂(t, ξ...) == equ_mod.B₂(t,ξ)
    @test equ_mod.B₃(t, ξ...) == equ_mod.B₃(t,ξ)

    @test equ_mod.b₁(t, ξ...) == equ_mod.b₁(t,ξ)
    @test equ_mod.b₂(t, ξ...) == equ_mod.b₂(t,ξ)
    @test equ_mod.b₃(t, ξ...) == equ_mod.b₃(t,ξ)

    @test equ_mod.E₁(t, ξ...) == equ_mod.E₁(t,ξ)
    @test equ_mod.E₂(t, ξ...) == equ_mod.E₂(t,ξ)
    @test equ_mod.E₃(t, ξ...) == equ_mod.E₃(t,ξ)

    @test equ_mod.A¹(t, ξ...) == equ_mod.A¹(t,ξ)
    @test equ_mod.A²(t, ξ...) == equ_mod.A²(t,ξ)
    @test equ_mod.A³(t, ξ...) == equ_mod.A³(t,ξ)

    @test equ_mod.B¹(t, ξ...) == equ_mod.B¹(t,ξ)
    @test equ_mod.B²(t, ξ...) == equ_mod.B²(t,ξ)
    @test equ_mod.B³(t, ξ...) == equ_mod.B³(t,ξ)

    @test equ_mod.b¹(t, ξ...) == equ_mod.b¹(t,ξ)
    @test equ_mod.b²(t, ξ...) == equ_mod.b²(t,ξ)
    @test equ_mod.b³(t, ξ...) == equ_mod.b³(t,ξ)

    @test equ_mod.E¹(t, ξ...) == equ_mod.E¹(t,ξ)
    @test equ_mod.E²(t, ξ...) == equ_mod.E²(t,ξ)
    @test equ_mod.E³(t, ξ...) == equ_mod.E³(t,ξ)

    @test equ_mod.dA₁dx₁(t, ξ...) == equ_mod.dA₁dx₁(t,ξ)
    @test equ_mod.dA₁dx₂(t, ξ...) == equ_mod.dA₁dx₂(t,ξ)
    @test equ_mod.dA₁dx₃(t, ξ...) == equ_mod.dA₁dx₃(t,ξ)

    @test equ_mod.dA₂dx₁(t, ξ...) == equ_mod.dA₂dx₁(t,ξ)
    @test equ_mod.dA₂dx₂(t, ξ...) == equ_mod.dA₂dx₂(t,ξ)
    @test equ_mod.dA₂dx₃(t, ξ...) == equ_mod.dA₂dx₃(t,ξ)

    @test equ_mod.dA₃dx₁(t, ξ...) == equ_mod.dA₃dx₁(t,ξ)
    @test equ_mod.dA₃dx₂(t, ξ...) == equ_mod.dA₃dx₂(t,ξ)
    @test equ_mod.dA₃dx₃(t, ξ...) == equ_mod.dA₃dx₃(t,ξ)
    
    @test equ_mod.ξ¹(t, ξ...) == equ_mod.ξ¹(t,ξ)
    @test equ_mod.ξ²(t, ξ...) == equ_mod.ξ²(t,ξ)
    @test equ_mod.ξ³(t, ξ...) == equ_mod.ξ³(t,ξ)
    
    @test equ_mod.ξ¹(t, ξ...) == equ_mod.ξ¹(t,ξ)
    @test equ_mod.ξ²(t, ξ...) == equ_mod.ξ²(t,ξ)
    @test equ_mod.ξ³(t, ξ...) == equ_mod.ξ³(t,ξ)
    
    # check internal consistency

    x = equ_mod.to_cartesian(t,ξ)

    @test equ_mod.from_cartesian(t, x) ≈ ξ  atol=1E-14

    g = equ_mod.g(t,ξ)
    ḡ = equ_mod.ḡ(t,ξ)
    DF = equ_mod.DF(t,ξ)
    DF̄ = equ_mod.DF̄(t,x)

    b = equ_mod.b(t,ξ)
    a⃗ = equ_mod.a⃗(t,ξ)
    b⃗ = equ_mod.b⃗(t,ξ)
    c⃗ = equ_mod.c⃗(t,ξ)


    @test equ_mod.J(t,ξ) ≈ sqrt(det(DF' * DF))  atol=1E-12
    @test ḡ         ≈ inv(g)          atol=1E-12
    @test DF' * DF  ≈ g               atol=1E-12
    @test DF  * DF̄  ≈ Array(I, 3, 3)  atol=1E-12
    @test DF̄  * DF̄' ≈ ḡ               atol=1E-12

    @test b⃗' * b     ≈ 1  atol=1E-14
    @test b' * ḡ * b ≈ 1  atol=1E-14
    @test a⃗' * g * a⃗ ≈ 1  atol=1E-14
    @test b⃗' * g * b⃗ ≈ 1  atol=1E-14
    @test c⃗' * g * c⃗ ≈ 1  atol=1E-14

    @test a⃗' * g * b⃗ ≈ 0  atol=1E-14
    @test a⃗' * g * c⃗ ≈ 0  atol=1E-14
    @test b⃗' * g * c⃗ ≈ 0  atol=1E-14

end


# equilibrium list (equilibrium, parameters, periodicity, module)
eqs = (
    (AxisymmetricTokamakCartesian,               (2., 3., 2.),   [0., 0., 0.]),
    (AxisymmetricTokamakCircular,                (2., 3., 2.),   [0., 2π, 2π]),
    (AxisymmetricTokamakCylindrical,             (2., 3., 2.),   [0., 0., 2π]),
    (AxisymmetricTokamakToroidalRegularization,  (2., 3., 2.),   [0., 2π, 2π]),
    (Singular,                                   (1.),           [0., 0., 0.]),
    (SymmetricQuadratic,                         (1.),           [0., 0., 0.]),
    (ThetaPinch,                                 (1.),           [0., 0., 0.]),
    (ABC,                                        (1., 0.5, 0.5), [0., 0., 0.]),
    (Solovev,           (6.2, 5.3, 0.32, 1.8, 0.45, -0.155),                [0., 0., 2π]),
    (SolovevXpoint,     (6.2, 5.3, 0.32, 1.8, 0.45, -0.155, 0.88, -0.60),   [0., 0., 2π]),
    (SolovevQuadratic,  (6.2, 5.3, 1., 1.),                                 [0., 0., 2π]),
)


# perturbation list (equilibrium, parameters, perturbation, parameters, module)
perts = (
    (SymmetricQuadratic,    (1.),   EzCosZ,    (2.)),
    (ThetaPinch,            (1.),   EzCosZ,    (2.)),
)


# testing parameters
t = 1.
x = [1.05, 0.5, 0.5]

# test equilibria
@testset "$(rpad(teststring(equ[1]),60))" for equ in eqs begin
        equ_obj = equ[1].init(equ[2]...)
        test_equilibrium(equ[1], t, x)
        @test periodicity(x, equ_obj) == equ[3]
    end
end
println()


# test perturbations
@testset "$(rpad(teststring(equ[1], equ[3]),60))" for equ in perts begin
        equ_obj = equ[1].init(equ[2]..., perturbation=equ[3].init(equ[4]...))
        test_equilibrium(equ[1], t, x)
    end
end
println()


# test correctness of some of the magnetic fields

function test_axisymmetric_tokamak_cartesian_equilibrium(equ_mod, t=0., x=[1.5, 0.0, 0.5])
    @test equ_mod.B¹(t,x) == equ_mod.B₁(t,x)
    @test equ_mod.B²(t,x) == equ_mod.B₂(t,x)
    @test equ_mod.B³(t,x) == equ_mod.B₃(t,x)

    @test equ_mod.B₁(t,x) == - equ_mod.B₀ / equ_mod.q₀ * (equ_mod.q₀ * equ_mod.R₀ * equ_mod.Y(t,x) + equ_mod.X(t,x) * equ_mod.Z(t,x) ) / equ_mod.R(t,x)^2
    @test equ_mod.B₂(t,x) == + equ_mod.B₀ / equ_mod.q₀ * (equ_mod.q₀ * equ_mod.R₀ * equ_mod.X(t,x) - equ_mod.Y(t,x) * equ_mod.Z(t,x) ) / equ_mod.R(t,x)^2
    @test equ_mod.B₃(t,x) == + equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)
end

function test_axisymmetric_tokamak_circular_equilibrium(equ_mod, t=0., x=[0.5, π/10, π/5])
    @test equ_mod.B¹(t,x) == 0
    @test equ_mod.B²(t,x) == - equ_mod.B₀ / equ_mod.q₀ / equ_mod.R(t,x)
    @test equ_mod.B³(t,x) ≈  - equ_mod.B₀ * equ_mod.R₀ / equ_mod.R(t,x)^2 atol=1E-14

    @test equ_mod.B₁(t,x) == 0
    @test equ_mod.B₂(t,x) == - equ_mod.B₀ / equ_mod.q₀ * equ_mod.r(t,x)^2 / equ_mod.R(t,x)
    @test equ_mod.B₃(t,x) ≈  - equ_mod.B₀ * equ_mod.R₀ atol=1E-14
end

function test_axisymmetric_tokamak_cylindrical_equilibrium(equ_mod, t=0., x=[1.5, 0.5, π/5])
    @test equ_mod.B¹(t,x) == + equ_mod.B₀ / equ_mod.q₀ * equ_mod.Z(t,x) / equ_mod.R(t,x)
    @test equ_mod.B²(t,x) == - equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)
    @test equ_mod.B³(t,x) == - equ_mod.B₀ * equ_mod.R₀ / equ_mod.R(t,x)^2

    @test equ_mod.B₁(t,x) == + equ_mod.B₀ / equ_mod.q₀ * equ_mod.Z(t,x) / equ_mod.R(t,x)
    @test equ_mod.B₂(t,x) == - equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)
    @test equ_mod.B₃(t,x) == - equ_mod.B₀ * equ_mod.R₀
end

function test_consistency_axisymmetric_tokamak_circular_equilibrium(equ_tor, equ_car, t=0., ξ=[0.5, π/10, π/5])
    x  = equ_tor.to_cartesian(t,ξ)
    DF = equ_tor.DF(t,ξ)

    B_tor = [equ_tor.B¹(t,ξ), equ_tor.B²(t,ξ), equ_tor.B³(t,ξ)]
    B_car = [equ_car.B¹(t,x), equ_car.B²(t,x), equ_car.B³(t,x)]

    B̂_tor = [equ_tor.B₁(t,ξ), equ_tor.B₂(t,ξ), equ_tor.B₃(t,ξ)]
    B̂_car = [equ_car.B₁(t,x), equ_car.B₂(t,x), equ_car.B₃(t,x)]

    @test B_tor' * B̂_tor ≈ B_car' * B̂_car atol=1E-12
end

function test_consistency_axisymmetric_tokamak_cylindrical_equilibrium(equ_cyl, equ_car, t=0., ξ=[1.5, 0.5, π/5])
    x = equ_cyl.to_cartesian(t,ξ)

    DF = equ_cyl.DF(t,ξ)

    B_cyl = [equ_cyl.B¹(t,ξ), equ_cyl.B²(t,ξ), equ_cyl.B³(t,ξ)]
    B_car = [equ_car.B¹(t,x), equ_car.B²(t,x), equ_car.B³(t,x)]

    B̂_cyl = [equ_cyl.B₁(t,ξ), equ_cyl.B₂(t,ξ), equ_cyl.B₃(t,ξ)]
    B̂_car = [equ_car.B₁(t,x), equ_car.B₂(t,x), equ_car.B₃(t,x)]

    @test B_cyl' * B̂_cyl ≈ B_car' * B̂_car atol=1E-12
end

function test_symmetric_quadratic_equilibrium(equ_mod, t=0., x=[1.0, 0.5, 0.5])
    @test equ_mod.B¹(t,x) == equ_mod.B₁(t,x)
    @test equ_mod.B²(t,x) == equ_mod.B₂(t,x)
    @test equ_mod.B³(t,x) == equ_mod.B₃(t,x)

    @test equ_mod.B₁(t,x) == 0
    @test equ_mod.B₂(t,x) == 0
    @test equ_mod.B₃(t,x) == equ_mod.B₀ * (1 + equ_mod.X(t,x)^2 + equ_mod.Y(t,x)^2)

    @test equ_mod.B(t,x)  == equ_mod.B₀ * (1 + equ_mod.X(t,x)^2 + equ_mod.Y(t,x)^2)

    @test equ_mod.b¹(t,x) == 0
    @test equ_mod.b²(t,x) == 0
    @test equ_mod.b³(t,x) == 1

    @test equ_mod.b₁(t,x) == 0
    @test equ_mod.b₂(t,x) == 0
    @test equ_mod.b₃(t,x) == 1
end

function test_theta_pinch_equilibrium(equ_mod, t=0., x=[1.0, 0.5, 0.5])
    @test equ_mod.B¹(t,x) == equ_mod.B₁(t,x)
    @test equ_mod.B²(t,x) == equ_mod.B₂(t,x)
    @test equ_mod.B³(t,x) == equ_mod.B₃(t,x)

    @test equ_mod.B₁(t,x) == 0
    @test equ_mod.B₂(t,x) == 0
    @test equ_mod.B₃(t,x) == equ_mod.B₀

    @test equ_mod.B(t,x)  == equ_mod.B₀

    @test equ_mod.b¹(t,x) == 0
    @test equ_mod.b²(t,x) == 0
    @test equ_mod.b³(t,x) == 1

    @test equ_mod.b₁(t,x) == 0
    @test equ_mod.b₂(t,x) == 0
    @test equ_mod.b₃(t,x) == 1
end

function test_abc_equilibrium(equ_mod, t=0., x=[1.0, 0.5, 0.5])
    @test equ_mod.B¹(t,x) == equ_mod.B₁(t,x)
    @test equ_mod.B²(t,x) == equ_mod.B₂(t,x)
    @test equ_mod.B³(t,x) == equ_mod.B₃(t,x)

    @test equ_mod.B₁(t,x) == equ_mod.A₁(t,x)
    @test equ_mod.B₂(t,x) == equ_mod.A₂(t,x)
    @test equ_mod.B₃(t,x) == equ_mod.A₃(t,x)

    @test equ_mod.B₁(t,x) == equ_mod.a₀ * sin(equ_mod.Z(t,x)) + equ_mod.c₀ * cos(equ_mod.Y(t,x))
    @test equ_mod.B₂(t,x) == equ_mod.b₀ * sin(equ_mod.X(t,x)) + equ_mod.a₀ * cos(equ_mod.Z(t,x))
    @test equ_mod.B₃(t,x) == equ_mod.c₀ * sin(equ_mod.Y(t,x)) + equ_mod.b₀ * cos(equ_mod.X(t,x))
end

@testset "$(rpad("Magnetic Fields",60))" begin
    test_axisymmetric_tokamak_cartesian_equilibrium(AxisymmetricTokamakCartesian)
    test_axisymmetric_tokamak_circular_equilibrium(AxisymmetricTokamakCircular)
    test_axisymmetric_tokamak_cylindrical_equilibrium(AxisymmetricTokamakCylindrical)
    test_symmetric_quadratic_equilibrium(SymmetricQuadratic)
    test_theta_pinch_equilibrium(ThetaPinch)
    test_abc_equilibrium(ABC)
end

@testset "$(rpad("Consistency",60))" begin
    test_consistency_axisymmetric_tokamak_circular_equilibrium(AxisymmetricTokamakCircular, AxisymmetricTokamakCartesian)
    test_consistency_axisymmetric_tokamak_cylindrical_equilibrium(AxisymmetricTokamakCylindrical, AxisymmetricTokamakCartesian)
end

println()
