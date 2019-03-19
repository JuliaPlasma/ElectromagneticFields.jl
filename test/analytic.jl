
function test_equilibrium(equ_mod, t=0., x=[2., 1., 4.])
    @test equ_mod.A₁(x...) == equ_mod.A₁(t,x)
    @test equ_mod.A₂(x...) == equ_mod.A₂(t,x)
    @test equ_mod.A₃(x...) == equ_mod.A₃(t,x)

    @test equ_mod.B₁(x...) == equ_mod.B₁(t,x)
    @test equ_mod.B₂(x...) == equ_mod.B₂(t,x)
    @test equ_mod.B₃(x...) == equ_mod.B₃(t,x)

    @test equ_mod.b₁(x...) == equ_mod.b₁(t,x)
    @test equ_mod.b₂(x...) == equ_mod.b₂(t,x)
    @test equ_mod.b₃(x...) == equ_mod.b₃(t,x)

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

module AxisymmetricTokamakCartesianEquilibrium end
module AxisymmetricTokamakCylindricalEquilibrium end
module AxisymmetricTokamakEquilibrium end
module SolovevEquilibrium end
module SolovevXpointEquilibrium end
module SolovevQuadraticEquilibrium end
module SymmetricQuadraticEquilibrium end
module ThetaPinchEquilibrium end
module ABCEquilibrium end

load_equilibrium(equ1, target_module=AxisymmetricTokamakCartesianEquilibrium)
load_equilibrium(equ2, target_module=AxisymmetricTokamakCylindricalEquilibrium)
load_equilibrium(equ3, target_module=AxisymmetricTokamakEquilibrium)
load_equilibrium(equ4, target_module=SolovevEquilibrium)
load_equilibrium(equ5, target_module=SolovevXpointEquilibrium)
load_equilibrium(equ6, target_module=SolovevQuadraticEquilibrium)
load_equilibrium(equ7, target_module=SymmetricQuadraticEquilibrium)
load_equilibrium(equ8, target_module=ThetaPinchEquilibrium)
load_equilibrium(equ9, target_module=ABCEquilibrium)

test_equilibrium(AxisymmetricTokamakCartesianEquilibrium)
test_equilibrium(AxisymmetricTokamakCylindricalEquilibrium)
test_equilibrium(AxisymmetricTokamakEquilibrium)
test_equilibrium(SolovevEquilibrium)
test_equilibrium(SolovevXpointEquilibrium)
test_equilibrium(SolovevQuadraticEquilibrium)
test_equilibrium(SymmetricQuadraticEquilibrium)
test_equilibrium(ThetaPinchEquilibrium)
test_equilibrium(ABCEquilibrium)
