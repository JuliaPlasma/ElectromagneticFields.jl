# Axisymmetric Tokamak (Cylindrical)

```@docs
AxisymmetricTokamakCylindrical
```

## Constructing the Field

```@example atcyl
using CairoMakie
using ElectromagneticFields

equ = AxisymmetricTokamakCylindrical.init()
```

## Plotting

In cylindrical coordinates the flux surfaces are the contours of the poloidal flux function
``A_\phi / R``:

```@example atcyl
plot_equilibrium(equ)
```

## Evaluating the Field

```@example atcyl
AxisymmetricTokamakCylindrical.@code()
nothing # hide
```

The generated functions take the time followed by the three coordinates ``(R, Z, \phi)``. Here
the absolute value of the magnetic field and the three components of the vector potential are
sampled on a poloidal grid that extends well beyond the plasma:

```@example atcyl
nr, nz = 100, 120

Rgrid = LinRange(0.25, 2.75, nr)
Zgrid = LinRange(-2.0, +2.0, nz)

sample(f) = [f(0.0, Rgrid[i], Zgrid[j], 0.0) for i in eachindex(Rgrid), j in eachindex(Zgrid)]

Bfield = sample(B)
A_R = sample(A₁)
A_Z = sample(A₂)
A_ϕ = sample(A₃)

extrema(Bfield)
```

The magnetic field is strongest on the inboard side, and it is ``A_\phi`` whose contours are the
flux surfaces. The other two components of the vector potential generate the toroidal field:
``A_R`` grows linearly in ``Z``, while ``A_Z`` depends on ``R`` alone, so its contours are
vertical lines.

```@example atcyl
fig = Figure(size = (800, 800))

panels = ((Bfield, L"|B| (R,Z)"), (A_ϕ, L"A_\phi (R,Z)"),
          (A_R, L"A_R (R,Z)"), (A_Z, L"A_Z (R,Z)"))

for (n, (vals, title)) in enumerate(panels)
    ax = Axis(fig[cld(n,2), mod1(n,2)];
        xlabel = L"R", ylabel = L"Z", title = title, aspect = DataAspect())
    contour!(ax, Rgrid, Zgrid, vals; levels = 20)
end

fig
```
