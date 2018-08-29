
function test_equilibrium(equ_mod, t=0., x=[2., 1., 4.])
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


equ1 = AxisymmetricTokamakCylindrical(1., 2., 2.)
equ2 = AxisymmetricTokamakToroidal(1., 2., 2.)
equ3 = Solovev(6.2, 5.3, 0.32, 1.8, 0.45, -0.155)
equ4 = SolovevXpoint(6.2, 5.3, 0.32, 1.8, 0.45, -0.155, 0.88, -0.60)

module AxisymmetricTokamakCylindricalEquilibrium end
module AxisymmetricTokamakEquilibrium end
module SolovevEquilibrium end
module SolovevXpointEquilibrium end

load_equilibrium(equ1, target_module=AxisymmetricTokamakCylindricalEquilibrium)
load_equilibrium(equ2, target_module=AxisymmetricTokamakEquilibrium)
load_equilibrium(equ3, target_module=SolovevEquilibrium)
load_equilibrium(equ4, target_module=SolovevXpointEquilibrium)

test_equilibrium(AxisymmetricTokamakCylindricalEquilibrium)
test_equilibrium(AxisymmetricTokamakEquilibrium)
test_equilibrium(SolovevEquilibrium)
test_equilibrium(SolovevXpointEquilibrium)
