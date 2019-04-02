
function test_equilibrium(equ_mod, t=0., x=[2., 1., 4.])
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


equ1 = AxisymmetricTokamakCartesian(1., 2., 2.)
equ2 = AxisymmetricTokamakCylindrical(1., 2., 2.)
equ3 = AxisymmetricTokamakToroidal(1., 2., 2.)
equ4 = Solovev(6.2, 5.3, 0.32, 1.8, 0.45, -0.155)
equ5 = SolovevXpoint(6.2, 5.3, 0.32, 1.8, 0.45, -0.155, 0.88, -0.60)
equ6 = SolovevQuadratic(6.2, 5.3, 1., 1.)
equ7 = SymmetricQuadratic(1.)
equ8 = ThetaPinch(1.)
equ9 = ABC(1., 0.5, 0.5)

pert1 = EzCosZ(2.)

module AxisymmetricTokamakCartesianEquilibrium end
module AxisymmetricTokamakCylindricalEquilibrium end
module AxisymmetricTokamakToroidalEquilibrium end
module SolovevEquilibrium end
module SolovevXpointEquilibrium end
module SolovevQuadraticEquilibrium end
module SymmetricQuadraticEquilibrium end
module SymmetricQuadraticEquilibriumEzCosZPerturbation end
module ThetaPinchEquilibrium end
module ThetaPinchEquilibriumEzCosZPerturbation end
module ABCEquilibrium end

load_equilibrium(equ1, target_module=AxisymmetricTokamakCartesianEquilibrium; output=1)
load_equilibrium(equ2, target_module=AxisymmetricTokamakCylindricalEquilibrium; output=1)
load_equilibrium(equ3, target_module=AxisymmetricTokamakToroidalEquilibrium; output=1)
load_equilibrium(equ4, target_module=SolovevEquilibrium; output=1)
load_equilibrium(equ5, target_module=SolovevXpointEquilibrium; output=1)
load_equilibrium(equ6, target_module=SolovevQuadraticEquilibrium; output=1)
load_equilibrium(equ7, target_module=SymmetricQuadraticEquilibrium; output=1)
load_equilibrium(equ8, target_module=ThetaPinchEquilibrium; output=1)
load_equilibrium(equ9, target_module=ABCEquilibrium; output=1)

load_equilibrium(equ7, pert1, target_module=SymmetricQuadraticEquilibriumEzCosZPerturbation; output=1)
load_equilibrium(equ8, pert1, target_module=ThetaPinchEquilibriumEzCosZPerturbation; output=1)

test_equilibrium(AxisymmetricTokamakCartesianEquilibrium)
test_equilibrium(AxisymmetricTokamakCylindricalEquilibrium)
test_equilibrium(AxisymmetricTokamakToroidalEquilibrium)
test_equilibrium(SolovevEquilibrium)
test_equilibrium(SolovevXpointEquilibrium)
test_equilibrium(SolovevQuadraticEquilibrium)
test_equilibrium(SymmetricQuadraticEquilibrium)
test_equilibrium(ThetaPinchEquilibrium)
test_equilibrium(ABCEquilibrium)

test_equilibrium(SymmetricQuadraticEquilibriumEzCosZPerturbation)
test_equilibrium(ThetaPinchEquilibriumEzCosZPerturbation)


x = [1., 1., 0.]

@test periodicity(x, equ1) == [0., 0., 0.]
@test periodicity(x, equ2) == [0., 0., 2π]
@test periodicity(x, equ3) == [0., 2π, 2π]
@test periodicity(x, equ4) == [0., 0., 2π]
@test periodicity(x, equ5) == [0., 0., 2π]
@test periodicity(x, equ6) == [0., 0., 2π]
@test periodicity(x, equ7) == [0., 0., 0.]
@test periodicity(x, equ8) == [0., 0., 0.]
@test periodicity(x, equ9) == [0., 0., 0.]


function test_axisymmetric_tokamak_cylindrical_equilibrium(equ_mod, t=0., x=[1.5, 0.5, π])
    @test equ_mod.B¹(t,x) ≈ - equ_mod.B₀ / equ_mod.q₀ * equ_mod.Z(t,x) / equ_mod.R(t,x)
    @test equ_mod.B²(t,x) ≈ equ_mod.B₀ / equ_mod.q₀ * (equ_mod.R(t,x) - equ_mod.R₀) / equ_mod.R(t,x)
    @test equ_mod.B³(t,x) ≈ - equ_mod.B₀ * equ_mod.R₀ / equ_mod.R(t,x)^2
end

function test_axisymmetric_tokamak_toroidal_equilibrium(equ_mod, t=0., x=[0.5, π/10, π])
    @test equ_mod.B¹(t,x) == zero(eltype(x))
    @test equ_mod.B²(t,x) ≈ equ_mod.B₀ / equ_mod.q₀ / equ_mod.R(t,x)
    @test equ_mod.B³(t,x) ≈ equ_mod.B₀ * equ_mod.R₀ / equ_mod.R(t,x)^2
end

println("Verifying magnetic field of axisymmetric tokamak equilibrium in cylindrical coordinates")
test_axisymmetric_tokamak_cylindrical_equilibrium(AxisymmetricTokamakCylindricalEquilibrium)

println("Verifying magnetic field of axisymmetric tokamak equilibrium in toroidal coordinates")
test_axisymmetric_tokamak_toroidal_equilibrium(AxisymmetricTokamakToroidalEquilibrium)
