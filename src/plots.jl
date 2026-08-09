
@doc raw"""
Plot an analytic equilibrium, typically as a contour plot of its vector potential.

```julia
plot_equilibrium(equ; size = ..., kwargs...)
```

Creates a new `Makie.Figure`, draws `equ` into it and returns the figure.

Plotting is provided by a package extension, so `Makie` (or one of its backends,
e.g. `CairoMakie` or `GLMakie`) needs to be loaded for these methods to exist:

```julia
using CairoMakie
using ElectromagneticFields

plot_equilibrium(Solovev.ITER())
```

Which quantity is shown depends on the equilibrium: for the tokamak and Solov'ev
equilibria it is the poloidal flux function ``A_\phi / R``, for the fields defined
in cartesian coordinates it is a panel per vector potential component, and for the
ABC field it is the absolute value of the magnetic field in three mid-planes.

All methods accept the resolution of the evaluation grid (`nx`, `ny`), the number of
contour `levels`, the plot ranges `xlims` and `ylims`, and the figure `size`.

See also [`plot_equilibrium!`](@ref).
"""
function plot_equilibrium end

@doc raw"""
Plot an analytic equilibrium into an existing figure.

```julia
plot_equilibrium!(position, equ; kwargs...)
```

Draws `equ` at `position`, which is any `Makie` grid position such as `fig[1,2]`,
and returns that position. This is what to use for composing several equilibria into
a single figure:

```julia
using CairoMakie
using ElectromagneticFields

fig = Figure(size = (1200, 400))
plot_equilibrium!(fig[1,1], Solovev.ITER())
plot_equilibrium!(fig[1,2], Solovev.NSTX(); xlims = (0.05, 2.3), ylims = (-2.25, +2.25))
fig
```

Like [`plot_equilibrium`](@ref), this requires `Makie` to be loaded.
"""
function plot_equilibrium! end
