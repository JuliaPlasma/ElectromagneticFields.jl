
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


# testing parameters
const t = 1.
const ξ = [1.05, 0.5, 0.5]


macro test_equilibrium(equilibrium_module, equilibrium_periodicity)
    module_name = string(equilibrium_module)
    module_test = Symbol(module_name * "Test")

    quote
        @eval module $module_test
            using Test
            using LinearAlgebra

            # import analytic field module
            import ElectromagneticFields
            import ElectromagneticFields.$equilibrium_module

            # inject code
            $equilibrium_module.@code

            function test_equilibrium(t, ξ)
                @test ξ¹(t, ξ...) == ξ¹(t,ξ)
                @test ξ²(t, ξ...) == ξ²(t,ξ)
                @test ξ³(t, ξ...) == ξ³(t,ξ)
                
                @test x¹(t, ξ...) == x¹(t,ξ)
                @test x²(t, ξ...) == x²(t,ξ)
                @test x³(t, ξ...) == x³(t,ξ)

                @test X(t, ξ...) == X(t,ξ)
                @test Y(t, ξ...) == Y(t,ξ)
                @test Z(t, ξ...) == Z(t,ξ)

                @test J(t, ξ...) == J(t,ξ)
                @test B(t, ξ...) == B(t,ξ)
                @test φ(t, ξ...) == φ(t,ξ)

                @test A₁(t, ξ...) == A₁(t,ξ)
                @test A₂(t, ξ...) == A₂(t,ξ)
                @test A₃(t, ξ...) == A₃(t,ξ)

                @test B₁(t, ξ...) == B₁(t,ξ)
                @test B₂(t, ξ...) == B₂(t,ξ)
                @test B₃(t, ξ...) == B₃(t,ξ)

                @test b₁(t, ξ...) == b₁(t,ξ)
                @test b₂(t, ξ...) == b₂(t,ξ)
                @test b₃(t, ξ...) == b₃(t,ξ)

                @test E₁(t, ξ...) == E₁(t,ξ)
                @test E₂(t, ξ...) == E₂(t,ξ)
                @test E₃(t, ξ...) == E₃(t,ξ)

                @test A¹(t, ξ...) == A¹(t,ξ)
                @test A²(t, ξ...) == A²(t,ξ)
                @test A³(t, ξ...) == A³(t,ξ)

                @test B¹(t, ξ...) == B¹(t,ξ)
                @test B²(t, ξ...) == B²(t,ξ)
                @test B³(t, ξ...) == B³(t,ξ)

                @test b¹(t, ξ...) == b¹(t,ξ)
                @test b²(t, ξ...) == b²(t,ξ)
                @test b³(t, ξ...) == b³(t,ξ)

                @test E¹(t, ξ...) == E¹(t,ξ)
                @test E²(t, ξ...) == E²(t,ξ)
                @test E³(t, ξ...) == E³(t,ξ)

                @test dA₁dx₁(t, ξ...) == dA₁dx₁(t,ξ)
                @test dA₁dx₂(t, ξ...) == dA₁dx₂(t,ξ)
                @test dA₁dx₃(t, ξ...) == dA₁dx₃(t,ξ)

                @test dA₂dx₁(t, ξ...) == dA₂dx₁(t,ξ)
                @test dA₂dx₂(t, ξ...) == dA₂dx₂(t,ξ)
                @test dA₂dx₃(t, ξ...) == dA₂dx₃(t,ξ)

                @test dA₃dx₁(t, ξ...) == dA₃dx₁(t,ξ)
                @test dA₃dx₂(t, ξ...) == dA₃dx₂(t,ξ)
                @test dA₃dx₃(t, ξ...) == dA₃dx₃(t,ξ)
                
                
                
                @test periodicity(t,ξ...) == $equilibrium_periodicity
                @test periodicity(t,ξ)    == $equilibrium_periodicity
                @test periodicity(ξ...)   == $equilibrium_periodicity
                @test periodicity(ξ)      == $equilibrium_periodicity

                # check internal consistency
                @test from_cartesian(t, to_cartesian(t,ξ)) ≈ ξ  atol=1E-14

                let g = g(t,ξ), ḡ = ḡ(t,ξ), DF = DF(t,ξ), DF̄ = DF̄(t,ξ),
                    a = a(t,ξ), b = b(t,ξ), c = c(t,ξ),
                    a⃗ = a⃗(t,ξ), b⃗ = b⃗(t,ξ), c⃗ = c⃗(t,ξ),
                    â = aₚ(t,ξ), b̂ = bₚ(t,ξ), ĉ = cₚ(t,ξ)

                    @test J(t,ξ) ≈ sqrt(det(DF' * DF))  atol=1E-12
                    @test ḡ         ≈ inv(g)          atol=1E-12
                    @test DF̄        ≈ inv(DF)         atol=1E-12
                    @test DF' * DF  ≈ g               atol=1E-12
                    @test DF  * DF̄  ≈ Array(I, 3, 3)  atol=1E-12
                    @test DF̄  * DF̄' ≈ ḡ               atol=1E-12


                    if $equilibrium_module != ElectromagneticFields.Singular
                        @test g * a⃗ ≈ a  atol=1E-14
                        @test g * b⃗ ≈ b  atol=1E-14
                        @test g * c⃗ ≈ c  atol=1E-14

                        @test ḡ * a ≈ a⃗  atol=1E-14
                        @test ḡ * b ≈ b⃗  atol=1E-14
                        @test ḡ * c ≈ c⃗  atol=1E-14

                        @test â ≈ DF * a⃗  atol=1E-14
                        @test b̂ ≈ DF * b⃗  atol=1E-14
                        @test ĉ ≈ DF * c⃗  atol=1E-14

                        @test â ≈ DF̄' * a  atol=1E-14
                        @test b̂ ≈ DF̄' * b  atol=1E-14
                        @test ĉ ≈ DF̄' * c  atol=1E-14

                        @test a⃗' * a ≈ 1  atol=1E-14
                        @test b⃗' * b ≈ 1  atol=1E-14
                        @test c⃗' * c ≈ 1  atol=1E-14

                        @test â' * â ≈ 1  atol=1E-14
                        @test b̂' * b̂ ≈ 1  atol=1E-14
                        @test ĉ' * ĉ ≈ 1  atol=1E-14

                        @test â' * b̂ ≈ 0  atol=1E-14
                        @test b̂' * ĉ ≈ 0  atol=1E-14
                        @test ĉ' * â ≈ 0  atol=1E-14

                        @test a' * ḡ * a ≈ 1  atol=1E-14
                        @test b' * ḡ * b ≈ 1  atol=1E-14
                        @test c' * ḡ * c ≈ 1  atol=1E-14

                        @test a' * ḡ * b ≈ 0  atol=1E-14
                        @test b' * ḡ * c ≈ 0  atol=1E-14
                        @test c' * ḡ * a ≈ 0  atol=1E-14

                        @test a⃗' * g * a⃗ ≈ 1  atol=1E-14
                        @test b⃗' * g * b⃗ ≈ 1  atol=1E-14
                        @test c⃗' * g * c⃗ ≈ 1  atol=1E-14

                        @test a⃗' * g * b⃗ ≈ 0  atol=1E-14
                        @test b⃗' * g * c⃗ ≈ 0  atol=1E-14
                        @test c⃗' * g * a⃗ ≈ 0  atol=1E-14
                    end
                end
            end
        end

        # run tests
        @testset "$(rpad($module_name,60))" begin
            import .$module_test
            $module_test.test_equilibrium(t, ξ)
        end
    end
end



# perturbation list (equilibrium, parameters, perturbation, parameters, module)
perts = (
    (SymmetricQuadratic,    EzCosZ,    (2.)),
    (ThetaPinch,            EzCosZ,    (2.)),
)




# test equilibria

@test_equilibrium ABC                                           [0., 0., 0.]
@test_equilibrium AxisymmetricTokamakCartesian                  [0., 0., 0.]
@test_equilibrium AxisymmetricTokamakCylindrical                [0., 0., 2π]
@test_equilibrium AxisymmetricTokamakToroidal                   [0., 2π, 2π]
@test_equilibrium AxisymmetricTokamakToroidalRegularization     [0., 2π, 2π]
@test_equilibrium PenningTrapUniform                            [0., 0., 0.]
@test_equilibrium PenningTrapBottle                             [0., 0., 0.]
@test_equilibrium PenningTrapAsymmetric                         [0., 0., 0.]
@test_equilibrium Singular                                      [0., 0., 0.]
@test_equilibrium SymmetricQuadratic                            [0., 0., 0.]
@test_equilibrium ThetaPinch                                    [0., 0., 0.]
# @test_equilibrium SolovevFRC                                    [0., 0., 2π]
@test_equilibrium SolovevITER                                   [0., 0., 2π]
@test_equilibrium SolovevITERwXpoint                            [0., 0., 2π]
@test_equilibrium SolovevNSTX                                   [0., 0., 2π]
@test_equilibrium SolovevNSTXwXpoint                            [0., 0., 2π]
@test_equilibrium SolovevNSTXwDoubleXpoint                      [0., 0., 2π]
@test_equilibrium SolovevSymmetric                              [0., 0., 0.]
println()


# test perturbations
# @testset "$(rpad(teststring(equ[1], equ[3]),60))" for equ in perts begin
#         equ_obj = equ[1].init(equ[2]..., perturbation=equ[3].init(equ[4]...))
#         # test_equilibrium(equ[1], t, x)
#     end
# end
# println()


# test correctness of some of the magnetic fields

function test_axisymmetric_tokamak_cartesian_equilibrium(equ_mod, t=0., x=[1.5, 0.0, 0.5])
    @test equ_mod.B¹(t,x) ≈ equ_mod.B₁(t,x)   atol=1E-16
    @test equ_mod.B²(t,x) ≈ equ_mod.B₂(t,x)   atol=1E-16
    @test equ_mod.B³(t,x) ≈ equ_mod.B₃(t,x)   atol=1E-16

    @test equ_mod.B₁(t,x) ≈ - equ_mod.B₀ / equ_mod.q₀ * (equ_mod.q₀ * equ_mod.R₀ * equ_mod.Y(t,x) + equ_mod.X(t,x) * equ_mod.Z(t,x) ) / equ_mod.R(t,x)^2   atol=1E-16
    @test equ_mod.B₂(t,x) ≈ + equ_mod.B₀ / equ_mod.q₀ * (equ_mod.q₀ * equ_mod.R₀ * equ_mod.X(t,x) - equ_mod.Y(t,x) * equ_mod.Z(t,x) ) / equ_mod.R(t,x)^2   atol=1E-16
    @test equ_mod.B₃(t,x) ≈ + equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)   atol=1E-16
end

function test_axisymmetric_tokamak_cylindrical_equilibrium(equ_mod, t=0., x=[1.5, 0.5, π/5])
    @test equ_mod.B¹(t,x) == + equ_mod.B₀ / equ_mod.q₀ * equ_mod.Z(t,x) / equ_mod.R(t,x)
    @test equ_mod.B²(t,x) == - equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)
    @test equ_mod.B³(t,x) == - equ_mod.B₀ * equ_mod.R₀ / equ_mod.R(t,x)^2

    @test equ_mod.B₁(t,x) == + equ_mod.B₀ / equ_mod.q₀ * equ_mod.Z(t,x) / equ_mod.R(t,x)
    @test equ_mod.B₂(t,x) == - equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)
    @test equ_mod.B₃(t,x) == - equ_mod.B₀ * equ_mod.R₀
end

function test_axisymmetric_tokamak_toroidal_equilibrium(equ_mod, t=0., x=[0.5, π/10, π/5])
    @test equ_mod.B¹(t,x) == 0
    @test equ_mod.B²(t,x) == - equ_mod.B₀ / equ_mod.q₀ / equ_mod.R(t,x)
    @test equ_mod.B³(t,x) ≈  - equ_mod.B₀ * equ_mod.R₀ / equ_mod.R(t,x)^2 atol=1E-14

    @test equ_mod.B₁(t,x) == 0
    @test equ_mod.B₂(t,x) == - equ_mod.B₀ / equ_mod.q₀ * equ_mod.r(t,x)^2 / equ_mod.R(t,x)
    @test equ_mod.B₃(t,x) ≈  - equ_mod.B₀ * equ_mod.R₀ atol=1E-14
end

function test_consistency_axisymmetric_tokamak_cylindrical_equilibrium(equ_cyl, equ_car, t=0., ξ=[1.5, 0.5, π/5])
    x  = equ_cyl.to_cartesian(t,ξ)
    DF = equ_cyl.DF(t,ξ)
    DF̄ = equ_cyl.DF̄(t,ξ)

    B_cyl = [equ_cyl.B¹(t,ξ), equ_cyl.B²(t,ξ), equ_cyl.B³(t,ξ)]
    B_car = [equ_car.B¹(t,x), equ_car.B²(t,x), equ_car.B³(t,x)]

    B̂_cyl = [equ_cyl.B₁(t,ξ), equ_cyl.B₂(t,ξ), equ_cyl.B₃(t,ξ)]
    B̂_car = [equ_car.B₁(t,x), equ_car.B₂(t,x), equ_car.B₃(t,x)]

    @test B_cyl' * B̂_cyl ≈ B_car' * B̂_car  atol=1E-12
    @test DF  * B_cyl ≈ - B_car  atol=1E-12
    @test DF̄' * B̂_cyl ≈ - B̂_car  atol=1E-12
end

function test_consistency_axisymmetric_tokamak_toroidal_equilibrium(equ_tor, equ_car, t=0., ξ=[0.5, π/10, π/5])
    x  = equ_tor.to_cartesian(t,ξ)
    DF = equ_tor.DF(t,ξ)
    DF̄ = equ_tor.DF̄(t,ξ)

    B_tor = [equ_tor.B¹(t,ξ), equ_tor.B²(t,ξ), equ_tor.B³(t,ξ)]
    B_car = [equ_car.B¹(t,x), equ_car.B²(t,x), equ_car.B³(t,x)]

    B̂_tor = [equ_tor.B₁(t,ξ), equ_tor.B₂(t,ξ), equ_tor.B₃(t,ξ)]
    B̂_car = [equ_car.B₁(t,x), equ_car.B₂(t,x), equ_car.B₃(t,x)]

    @test B_tor' * B̂_tor ≈ B_car' * B̂_car  atol=1E-12
    @test DF  * B_tor ≈ - B_car  atol=1E-12
    @test DF̄' * B̂_tor ≈ - B̂_car  atol=1E-12
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
    test_axisymmetric_tokamak_cartesian_equilibrium(AxisymmetricTokamakCartesianTest)
    test_axisymmetric_tokamak_cylindrical_equilibrium(AxisymmetricTokamakCylindricalTest)
    test_axisymmetric_tokamak_toroidal_equilibrium(AxisymmetricTokamakToroidalTest)
    test_symmetric_quadratic_equilibrium(SymmetricQuadraticTest)
    test_theta_pinch_equilibrium(ThetaPinchTest)
    test_abc_equilibrium(ABCTest)
end

@testset "$(rpad("Consistency",60))" begin
    test_consistency_axisymmetric_tokamak_cylindrical_equilibrium(AxisymmetricTokamakCylindricalTest, AxisymmetricTokamakCartesianTest)
    test_consistency_axisymmetric_tokamak_toroidal_equilibrium(AxisymmetricTokamakToroidalTest, AxisymmetricTokamakCartesianTest)
end

println()
