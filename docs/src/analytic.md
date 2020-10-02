
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

