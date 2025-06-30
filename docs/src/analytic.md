
# Analytic Fields


## Arnold-Beltrami-Childress (ABC) Field

```@docs
ABC
```
```@eval
using Plots
using ElectromagneticFields

eq_abc = ElectromagneticFields.ABC.init()
plot(eq_abc)
savefig("abc.png")

nothing
```

Absolut value of the magnetic field:

![](abc.png)



## Axisymmetric Tokamak Equilibrium in Cartesian Coordinates

```@docs
AxisymmetricTokamakCartesian
```
```@eval
using Plots
using ElectromagneticFields

eq_car = ElectromagneticFields.AxisymmetricTokamakCartesian.init()
plot(eq_car)
savefig("axisymmetric_tokamak_cartesian.png")

nothing
```

Vector potential in y direction:

![](axisymmetric_tokamak_cartesian.png)



## Axisymmetric Tokamak Equilibrium in Cylindrical Coordinates

```@docs
AxisymmetricTokamakCylindrical
```
```@eval
using Plots
using ElectromagneticFields

eq_cyl = ElectromagneticFields.AxisymmetricTokamakCylindrical.init()
plot(eq_cyl)
savefig("axisymmetric_tokamak_cylindrical.png")

nothing
```

Vector potential in y direction:

![](axisymmetric_tokamak_cylindrical.png)



## Axisymmetric Tokamak Equilibrium in Toroidal Coordinates

```@docs
AxisymmetricTokamakToroidal
```
```@eval
using Plots
using ElectromagneticFields

eq_cir = ElectromagneticFields.AxisymmetricTokamakToroidal.init()
plot(eq_cir)
savefig("axisymmetric_tokamak_toroidal.png")

nothing
```

Vector potential in y direction:

![](axisymmetric_tokamak_toroidal.png)



## Dipole
```@docs
Dipole
```
```@eval
using Plots
using ElectromagneticFields

eq_dipole = ElectromagneticFields.Dipole.init()
plot(eq_dipole)
savefig("dipole.png")

nothing
```

Vector potential components:

![](dipole.png)



## Penning Trap with Uniform Magnetic Field

```@docs
PenningTrapUniform
```



## Penning Trap with Magnetic Bottle

```@docs
PenningTrapBottle
```



## Penning Trap with Asymmetric Magnetic Field

```@docs
PenningTrapAsymmetric
```



## Quadratic Potentials
```@docs
QuadraticPotentials
```
```@eval
using Plots
using ElectromagneticFields

eq_quadratic = ElectromagneticFields.QuadraticPotentials.init()
plot(eq_quadratic)
savefig("quadratic-potentials.png")

nothing
```

Vector potential components:

![](quadratic-potentials.png)



## Singular Magnetic Field

```@docs
Singular
```
```@eval
using Plots
using ElectromagneticFields

eq_sng = ElectromagneticFields.Singular.init()
plot(eq_sng)
savefig("singular.png")

nothing
```

Vector potential and magnetic field components:

![](singular.png)



## Symmetric Magnetic Field

```@docs
SymmetricQuadratic
```
```@eval
using Plots
using ElectromagneticFields

eq_sym = ElectromagneticFields.SymmetricQuadratic.init()
plot(eq_sym)
savefig("symmetric.png")

nothing
```

Vector potential and magnetic field components:

![](symmetric.png)



## Symmetric Solov'ev Equilibrium

```@docs
SolovevSymmetric
```

```@eval
using Plots
using ElectromagneticFields

eq_sol = ElectromagneticFields.SolovevSymmetric.init(0.0, 1.0, 2.0, 0.5)
plot(eq_sol)
savefig("solovev_symmetric.png")

nothing
```

Vector potential in z direction:

![](solovev_symmetric.png)



## Solov'ev Equilibrium (Antoine Cerfon)

### Solov'ev Equilibrium up/down symmetric

```@docs
ElectromagneticFields.Solovev.SolovevEquilibrium
```

```@eval
using Plots
using ElectromagneticFields

eq_sol_iter = ElectromagneticFields.Solovev.ITER()
plot_iter = plot(eq_sol_iter, title="ITER", xlims=(0.6,1.4))

eq_sol_nstx = ElectromagneticFields.Solovev.NSTX()
plot_nstx = plot(eq_sol_nstx, title="NSTX", xlims=(0.05, 2.3), ylims=(-2.25, +2.25), levels=50)

eq_sol_frc = ElectromagneticFields.Solovev.FRC()
plot_frc = plot(eq_sol_frc, title="FRC", xlims=(0,2), ylims=(-10,+10), levels=25, aspect_ratio=1//5)

plot(size=(800,400), layout=(1,3), plot_iter, plot_nstx, plot_frc)
savefig("solovev.png")

nothing
```

Vector potential in ITER, NSTX and a field reversed configuration:

![](solovev.png)



### Solov'ev Equilibrium with X-Point

```@docs
ElectromagneticFields.Solovev.SolovevXpointEquilibrium
```

```@eval
using Plots
using ElectromagneticFields

eq_sol_iter = ElectromagneticFields.Solovev.ITER(xpoint=true)
plot_iter = plot(eq_sol_iter, title="ITER", xlims=(0.6,1.4))

eq_sol_nstx = ElectromagneticFields.Solovev.NSTX(xpoint=true)
plot_nstx = plot(eq_sol_nstx, title="NSTX", xlims=(0.05, 2.3), ylims=(-2.25, +2.25), levels=50)

eq_sol_dblx = ElectromagneticFields.Solovev.NSTXdoubleX()
plot_dblx = plot(eq_sol_dblx, title="NSTX", xlims=(0.05, 2.3), ylims=(-2.25, +2.25), levels=50)

plot(size=(800,400), layout=(1,3), plot_iter, plot_nstx, plot_dblx)
savefig("solovev_xpoint.png")

nothing
```

Vector potential in ITER and NSTX:

![](solovev_xpoint.png)



## Theta Pinch

```@docs
ThetaPinch
```
```@eval
using Plots
using ElectromagneticFields

eq_thp = ElectromagneticFields.ThetaPinch.init()
plot(eq_thp)
savefig("theta_pinch.png")

nothing
```

Vector potential components:

![](theta_pinch.png)
