
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

![](axisymmetric_tokamak_cylindrical.png)


## Axisymmetric Tokamak Equilibrium in Circular Coordinates

```@docs
AxisymmetricTokamakCircular
```
```@eval
using Plots
using ElectromagneticFields

eq_cir = ElectromagneticFields.AxisymmetricTokamakCircular.init()
plot(eq_cir)
savefig("axisymmetric_tokamak_circular.png")

nothing
```

![](axisymmetric_tokamak_circular.png)


