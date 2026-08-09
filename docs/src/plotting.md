# Plotting

Every analytic equilibrium can be plotted. The plotting routines are provided by a package
extension, so they become available as soon as [Makie](https://docs.makie.org) or one of its
backends is loaded. For documentation and other static output CairoMakie is the natural choice:

```@example plotting
using CairoMakie
using ElectromagneticFields

plot_equilibrium(Solovev.ITER())
```

What is shown depends on the field. For the tokamak and Solov'ev equilibria it is the poloidal
flux function ``A_\phi / R``, whose contours are the flux surfaces, with the plasma boundary drawn
on top in red. For the fields defined in cartesian coordinates it is one panel per component of
the vector potential, and for the ABC field the absolute value of the magnetic field in the three
mid-planes.


## Adjusting the Plot

All methods accept the resolution of the evaluation grid (`nx`, `ny`), the number of contour
`levels`, the plot ranges `xlims` and `ylims`, the figure `size`, and a `title`. Anything not
recognised is forwarded to Makie's `contour!`, so e.g. `colormap` and `linewidth` work too:

```@example plotting
plot_equilibrium(Solovev.NSTX();
    xlims = (0.05, 2.3),
    ylims = (-2.25, +2.25),
    levels = 25,
    size = (350, 500),
    colormap = :plasma,
)
```

Contour panels can also be given a colorbar, which is off by default:

```@example plotting
plot_equilibrium(ThetaPinch.init(); colorbar = true, size = (900, 400))
```


## Composing Figures

`plot_equilibrium` creates a figure of its own. To draw into an existing one — to compare
several configurations side by side, or to combine an equilibrium with other plots — use
`plot_equilibrium!`, which takes any Makie grid position as its first argument:

```@example plotting
fig = Figure(size = (900, 400))

plot_equilibrium!(fig[1,1], Solovev.ITER();
    title = "ITER", xlims = (0.6, 1.4))
plot_equilibrium!(fig[1,2], Solovev.NSTX();
    title = "NSTX", xlims = (0.05, 2.3), ylims = (-2.25, +2.25))
plot_equilibrium!(fig[1,3], Solovev.FRC();
    title = "FRC", xlims = (0.0, 2.0), ylims = (-10.0, +10.0), levels = 25,
    aspect = AxisAspect(0.5))

fig
```

Note the `aspect` keyword on the last panel. By default the axes use `DataAspect()`, so that one
unit in ``R`` has the same length as one unit in ``Z``. That is the right choice almost
everywhere, but a field reversed configuration is so elongated that it would leave the panel a
thin sliver, hence the explicit aspect ratio.


## Plotting by Hand

The routines above cover the flux surfaces and the vector potential, which is what one usually
wants to look at. For everything else, evaluate the field yourself and plot the result — the
generated functions return plain numbers and the grids are plain arrays, so there is nothing
special about it:

```@example plotting
AxisymmetricTokamakCylindrical.@code()
nothing # hide
```

```@example plotting
Rgrid = LinRange(0.5, 1.5, 100)
Zgrid = LinRange(-0.5, 0.5, 120)

Bfield = [B(0.0, Rgrid[i], Zgrid[j], 0.0) for i in eachindex(Rgrid), j in eachindex(Zgrid)]

fig = Figure(size = (500, 400))
ax = Axis(fig[1,1]; xlabel = L"R", ylabel = L"Z", title = L"|B| (R,Z)", aspect = DataAspect())
cf = contourf!(ax, Rgrid, Zgrid, Bfield; levels = 20)
Colorbar(fig[1,2], cf)

fig
```

The index order matters: `Bfield[i,j]` has to hold the value at `(Rgrid[i], Zgrid[j])`, which is
what the comprehension above produces and what Makie expects.


## Reference

```@docs
plot_equilibrium
plot_equilibrium!
```
