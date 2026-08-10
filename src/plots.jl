
@doc raw"""
Plot an analytic equilibrium, typically as a contour plot of its vector potential.

```julia
plot_equilibrium(equ; size = ..., figure = NamedTuple(), kwargs...)
```

Creates a new `Makie.Figure`, draws `equ` into it and returns the figure.

Plotting is provided by a package extension, so `Makie` (or one of its backends,
e.g. `CairoMakie` or `GLMakie`) needs to be loaded for these methods to exist:

```julia
using CairoMakie
using ElectromagneticFields

plot_equilibrium(Solovev.ITER())
```

Which quantity is shown depends on the equilibrium. For the tokamak and Solov'ev
equilibria it is the poloidal flux function ``\psi``, whose contours are the flux
surfaces; that is the third covariant component `A₃` of the vector potential —
``A_\phi`` in the toroidal charts, ``R \, A_y`` in the cartesian tokamak. For the
remaining fields defined in cartesian coordinates it is a panel per vector potential
component, and for the ABC field the absolute value of the magnetic field in three
mid-planes.

The figure `size` defaults to one chosen per equilibrium, and everything in `figure`
is passed on to `Makie.Figure`. Both can carry a size; the precedence is the
per-equilibrium default first, then `figure`, then an explicit `size`. All remaining
keyword arguments are forwarded to [`plot_equilibrium!`](@ref).

See also [`plot_equilibrium!`](@ref).
"""
function plot_equilibrium end

@doc raw"""
Plot an analytic equilibrium into an existing figure.

```julia
plot_equilibrium!(position, equ; kwargs...)
```

Draws `equ` at `position`, which is any `Makie` grid position such as `fig[1,2]`, and
returns the `Axis` it created, or the vector of axes for the equilibria that draw more
than one panel. This is what to use for composing several equilibria into a single
figure:

```julia
using CairoMakie
using ElectromagneticFields

fig = Figure(size = (1200, 400))
plot_equilibrium!(fig[1,1], Solovev.ITER())
plot_equilibrium!(fig[1,2], Solovev.NSTX(); xlims = (0.05, 2.3), ylims = (-2.25, +2.25))
fig
```

Like [`plot_equilibrium`](@ref), this requires `Makie` to be loaded.

# Keyword Arguments

Common to every method:

  - `levels`: number of contour levels, or the levels themselves
  - `title`, `xlabel`, `ylabel`: axis labels, each panel providing its own default
  - `aspect`: axis aspect ratio, `DataAspect()` by default
  - `colorbar`: draw a colorbar next to each panel, `false` by default

Anything not recognised is forwarded to `Makie.contour!`, so e.g. `colormap` and
`linewidth` work as well.

The remaining keywords depend on the field. Most methods take the resolution of the
evaluation grid as `nx` and `ny` and the plot range as `xlims` and `ylims`. The
exceptions are

  - `ABCEquilibrium`, which is sampled on a cubic grid of `nx` points per direction
    covering ``[0, 2\pi]``, and takes `ni` for the index of the mid-plane, by default
    the grid point closest to ``\pi`` (exactly ``\pi`` for odd `nx`);
  - `SolovevEquilibrium`, which takes `boundary` to switch off the plasma boundary
    drawn in red on top of the flux surfaces, resolved with `nτ` points. It is the
    only equilibrium that draws one.

For the Solov'ev equilibria, an integer `levels` counts levels spaced uniformly in
``\psi`` and anchored to the flux on the magnetic axis, with one level on the plasma
boundary and a quarter of them inside it. Spreading them over the sampled range
instead would spend nearly all of them on the far field, where ``\psi`` grows without
bound. Passing the levels themselves bypasses this, as everywhere else.

Not every equilibrium has a plotting method: the three Penning traps have none, and
report so in an `ArgumentError`. See the Plotting page of the documentation for how
to sample and plot such a field directly.
"""
function plot_equilibrium! end
