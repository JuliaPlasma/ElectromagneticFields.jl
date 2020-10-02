
# Analytic Fields


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


