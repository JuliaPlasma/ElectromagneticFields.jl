
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


function test_equilibrium(equ_mod, t, x)
    @test equ_mod.X(x...) == equ_mod.X(t,x)
    @test equ_mod.Y(x...) == equ_mod.Y(t,x)
    @test equ_mod.Z(x...) == equ_mod.Z(t,x)

    @test equ_mod.J(x...) == equ_mod.J(t,x)
    @test equ_mod.B(x...) == equ_mod.B(t,x)
    @test equ_mod.φ(x...) == equ_mod.φ(t,x)

    @test equ_mod.A₁(x...) == equ_mod.A₁(t,x)
    @test equ_mod.A₂(x...) == equ_mod.A₂(t,x)
    @test equ_mod.A₃(x...) == equ_mod.A₃(t,x)

    @test equ_mod.B₁(x...) == equ_mod.B₁(t,x)
    @test equ_mod.B₂(x...) == equ_mod.B₂(t,x)
    @test equ_mod.B₃(x...) == equ_mod.B₃(t,x)

    @test equ_mod.b₁(x...) == equ_mod.b₁(t,x)
    @test equ_mod.b₂(x...) == equ_mod.b₂(t,x)
    @test equ_mod.b₃(x...) == equ_mod.b₃(t,x)

    @test equ_mod.E₁(x...) == equ_mod.E₁(t,x)
    @test equ_mod.E₂(x...) == equ_mod.E₂(t,x)
    @test equ_mod.E₃(x...) == equ_mod.E₃(t,x)

    @test equ_mod.A¹(x...) == equ_mod.A¹(t,x)
    @test equ_mod.A²(x...) == equ_mod.A²(t,x)
    @test equ_mod.A³(x...) == equ_mod.A³(t,x)

    @test equ_mod.B¹(x...) == equ_mod.B¹(t,x)
    @test equ_mod.B²(x...) == equ_mod.B²(t,x)
    @test equ_mod.B³(x...) == equ_mod.B³(t,x)

    @test equ_mod.b¹(x...) == equ_mod.b¹(t,x)
    @test equ_mod.b²(x...) == equ_mod.b²(t,x)
    @test equ_mod.b³(x...) == equ_mod.b³(t,x)

    @test equ_mod.E¹(x...) == equ_mod.E¹(t,x)
    @test equ_mod.E²(x...) == equ_mod.E²(t,x)
    @test equ_mod.E³(x...) == equ_mod.E³(t,x)

    @test equ_mod.dA₁dx₁(x...) == equ_mod.dA₁dx₁(t,x)
    @test equ_mod.dA₁dx₂(x...) == equ_mod.dA₁dx₂(t,x)
    @test equ_mod.dA₁dx₃(x...) == equ_mod.dA₁dx₃(t,x)

    @test equ_mod.dA₂dx₁(x...) == equ_mod.dA₂dx₁(t,x)
    @test equ_mod.dA₂dx₂(x...) == equ_mod.dA₂dx₂(t,x)
    @test equ_mod.dA₂dx₃(x...) == equ_mod.dA₂dx₃(t,x)

    @test equ_mod.dA₃dx₁(x...) == equ_mod.dA₃dx₁(t,x)
    @test equ_mod.dA₃dx₂(x...) == equ_mod.dA₃dx₂(t,x)
    @test equ_mod.dA₃dx₃(x...) == equ_mod.dA₃dx₃(t,x)
end


# modules that will hold the generated equilibrium code
module AxisymmetricTokamakCartesianEquilibrium end
module AxisymmetricTokamakCylindricalEquilibrium end
module AxisymmetricTokamakToroidalEquilibrium end
module SymmetricQuadraticEquilibrium end
module ThetaPinchEquilibrium end
module ABCEquilibrium end
module SolovevEquilibrium end
module SolovevXpointEquilibrium end
module SolovevQuadraticEquilibrium end

module SymmetricQuadraticEquilibriumEzCosZPerturbation end
module ThetaPinchEquilibriumEzCosZPerturbation end


# equilibrium list (equilibrium, parameters, periodicity, module)
eqs = (
    (AxisymmetricTokamakCartesian,      (2., 3., 2.),   [0., 0., 0.],   AxisymmetricTokamakCartesianEquilibrium),
    (AxisymmetricTokamakCylindrical,    (2., 3., 2.),   [0., 0., 2π],   AxisymmetricTokamakCylindricalEquilibrium),
    (AxisymmetricTokamakToroidal,       (2., 3., 2.),   [0., 2π, 2π],   AxisymmetricTokamakToroidalEquilibrium),
    (SymmetricQuadratic,                (1.),           [0., 0., 0.],   SymmetricQuadraticEquilibrium),
    (ThetaPinch,                        (1.),           [0., 0., 0.],   ThetaPinchEquilibrium),
    (ABC,                               (1., 0.5, 0.5), [0., 0., 0.],   ABCEquilibrium),
    # (Solovev,           (6.2, 5.3, 0.32, 1.8, 0.45, -0.155),                [0., 0., 2π],   SolovevEquilibrium),
    # (SolovevXpoint,     (6.2, 5.3, 0.32, 1.8, 0.45, -0.155, 0.88, -0.60),   [0., 0., 2π],   SolovevXpointEquilibrium),
    # (SolovevQuadratic,  (6.2, 5.3, 1., 1.),                                 [0., 0., 2π],   SolovevQuadraticEquilibrium),
)


# perturbation list (equilibrium, parameters, perturbation, parameters, module)
perts = (
    (SymmetricQuadratic,    (1.),   EzCosZ,    (2.),   SymmetricQuadraticEquilibriumEzCosZPerturbation),
    (ThetaPinch,            (1.),   EzCosZ,    (2.),   ThetaPinchEquilibriumEzCosZPerturbation),
)


# testing parameters
t = 1.
x = [1., 1., 1.]

# test equilibria
@testset "$(rpad(teststring(equ[1]),60))" for equ in eqs begin
        equ_obj = equ[1](equ[2]...)
        load_equilibrium(equ_obj, target_module=equ[4])
        test_equilibrium(equ[4], t, x)
        @test periodicity(x, equ_obj) == equ[3]
    end
end
println()


# test perturbations
@testset "$(rpad(teststring(equ[1], equ[3]),60))" for equ in perts begin
        equ_obj = equ[1](equ[2]...)
        prt_obj = equ[3](equ[4]...)
        load_equilibrium(equ_obj, prt_obj, target_module=equ[5])
        test_equilibrium(equ[5], t, x)
    end
end
println()


# test correctness of some of the magnetic fields

function test_axisymmetric_tokamak_cylindrical_equilibrium(equ_mod, t=0., x=[1.5, 0.5, π])
    @test equ_mod.B¹(t,x) == - equ_mod.B₀ / equ_mod.q₀ * equ_mod.Z(t,x) / equ_mod.R(t,x)
    @test equ_mod.B²(t,x) == + equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)
    @test equ_mod.B³(t,x) == - equ_mod.B₀ * equ_mod.R₀ / equ_mod.R(t,x)^2

    @test equ_mod.B₁(t,x) == - equ_mod.B₀ / equ_mod.q₀ * equ_mod.Z(t,x) / equ_mod.R(t,x)
    @test equ_mod.B₂(t,x) == + equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)
    @test equ_mod.B₃(t,x) == - equ_mod.B₀ * equ_mod.R₀
end

function test_axisymmetric_tokamak_toroidal_equilibrium(equ_mod, t=0., x=[0.5, π/10, π])
    @test equ_mod.B¹(t,x) == zero(eltype(x))
    @test equ_mod.B²(t,x) == equ_mod.B₀ / equ_mod.q₀ / equ_mod.R(t,x)
    @test equ_mod.B³(t,x) == equ_mod.B₀ * equ_mod.R₀ / equ_mod.R(t,x)^2

    @test equ_mod.B₁(t,x) == zero(eltype(x))
    @test equ_mod.B₂(t,x) == equ_mod.B₀ / equ_mod.q₀ * equ_mod.r(t,x)^2 / equ_mod.R(t,x)
    @test equ_mod.B₃(t,x) == equ_mod.B₀ * equ_mod.R₀
end

function test_theta_pinch_equilibrium(equ_mod, t=0., x=[1.5, 0.5, π])
    @test equ_mod.B¹(t,x) == 0
    @test equ_mod.B²(t,x) == 0
    @test equ_mod.B³(t,x) == equ_mod.B₀

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


@testset "$(rpad("Magnetic Fields",60))" begin
    test_axisymmetric_tokamak_cylindrical_equilibrium(AxisymmetricTokamakCylindricalEquilibrium)
    test_axisymmetric_tokamak_toroidal_equilibrium(AxisymmetricTokamakToroidalEquilibrium)
    test_theta_pinch_equilibrium(ThetaPinchEquilibrium)
end
println()
