

t = AxisymmetricTokamak(1., 2., 2.)
t = AxisymmetricTokamakToroidal(1., 2., 2.)
t = Solovev(6.2, 5.3, 0.32, 1.8, 0.45, -0.155)
t = SolovevXpoint(6.2, 5.3, 0.32, 1.8, 0.45, -0.155, 0.88, -0.60)

load_equilibrium(t)

@test dA₁dx₂(2., 1., 4.) == dA₁dx₂(0., [2., 1., 4.])
