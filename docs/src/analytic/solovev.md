# Solov'ev Equilibrium

The Solov'ev equilibria are exact solutions of the Grad-Shafranov equation, parametrised by the
inverse aspect ratio ``\epsilon``, the elongation ``\kappa`` and the triangularity ``\delta`` of
the plasma boundary. Between them these three shape parameters cover configurations as different
as a conventional tokamak, a spherical torus and a field reversed configuration.

The coordinates are ``(R/R_0, Z/R_0, \phi)``, i.e. normalised to the major radius.


## Up/Down Symmetric Equilibrium

```@docs
SolovevEquilibrium
```

### Constructing the Field

The parameters can be given explicitly,

```@example solovev
using CairoMakie
using ElectromagneticFields

equ = SolovevEquilibrium(6.2, 5.3, 0.32, 1.7, 0.33, -0.155)
```

or one of the named configurations can be used. `SolovevEquilibriumITER`,
`SolovevEquilibriumNSTX` and `SolovevEquilibriumFRC` are provided:

```@example solovev
equ = SolovevEquilibriumITER()
```

### Plotting

The contours are the flux surfaces, and the plasma boundary is drawn on top in red:

```@example solovev
plot_equilibrium(equ)
```

The three named configurations side by side. Note how much of the aspect ratio, elongation and
triangularity is visible directly in the shape of the boundary:

```@example solovev
fig = Figure(size = (900, 400))

plot_equilibrium!(fig[1,1], SolovevEquilibriumITER();
    title = "ITER", xlims = (0.6, 1.4))
plot_equilibrium!(fig[1,2], SolovevEquilibriumNSTX();
    title = "NSTX", xlims = (0.05, 2.3), ylims = (-2.25, +2.25))
plot_equilibrium!(fig[1,3], SolovevEquilibriumFRC();
    title = "FRC", xlims = (0.0, 2.0), ylims = (-10.0, +10.0),
    aspect = AxisAspect(0.5))

fig
```


## Equilibrium with X-Point

```@docs
SolovevXpointEquilibrium
```

### Constructing the Field

The `Xpoint` configurations place an X-point below the plasma, turning the outermost closed flux
surface into a separatrix:

```@example solovev
equ_xpoint = SolovevXpointEquilibriumITER()
```

[`SolovevDoubleXpointEquilibriumNSTX`](@ref) gives a configuration with X-points above *and*
below the plasma.

### Plotting

```@example solovev
fig = Figure(size = (900, 400))

plot_equilibrium!(fig[1,1], SolovevXpointEquilibriumITER();
    title = "ITER", xlims = (0.6, 1.4))
plot_equilibrium!(fig[1,2], SolovevXpointEquilibriumNSTX();
    title = "NSTX", xlims = (0.05, 2.3), ylims = (-2.25, +2.25))
plot_equilibrium!(fig[1,3], SolovevDoubleXpointEquilibriumNSTX();
    title = "NSTX (double X-point)", xlims = (0.05, 2.3), ylims = (-2.25, +2.25))

fig
```

The X-point shows up as the saddle where the innermost open surface pinches off — for the double
X-point configuration there is one at the top and one at the bottom.


## Evaluating the Field

```@example solovev
field = FieldFunctions(SolovevEquilibriumITER())
nothing # hide
```

The coordinates are normalised to ``R_0``, so a physical grid has to be divided by it before the
generated functions are called:

```@example solovev
R₀ = parameters(field).R₀

nr, nz = 100, 120

Rgrid = LinRange(3.0, 9.0, nr)
Zgrid = LinRange(-5.0, +5.0, nz)

sample(f) = [f(0.0, Rgrid[i] / R₀, Zgrid[j] / R₀, 0.0)
             for i in eachindex(Rgrid), j in eachindex(Zgrid)]

Bfield = sample((t, ξ...) -> B(field, t, ξ...))
A_R = sample((t, ξ...) -> A♭(field, t, ξ...)[1])
A_Z = sample((t, ξ...) -> A♭(field, t, ξ...)[2])
A_ϕ = sample((t, ξ...) -> A♭(field, t, ξ...)[3])

extrema(Bfield)
```

The plasma boundary is the flux surface parametrised by the shape parameters, which the field
carries as data:

```@example solovev
ϵ, κ, δ = parameters(field).ϵ, parameters(field).κ, parameters(field).δ

τ = LinRange(0, 2π, 200)

boundary_R = R₀ .* (1 .+ ϵ .* cos.(τ .+ asin(δ) .* sin.(τ)))
boundary_Z = R₀ .* ϵ .* κ .* sin.(τ)

nothing # hide
```

Putting the two together, the field strength and the three components of the vector potential
with the boundary drawn on top. The toroidal field dominates by far, so the contours of ``|B|``
are almost the lines of constant ``R``, and it is ``A_\phi`` that carries the flux surfaces:

```@example solovev
fig = Figure(size = (800, 900))

panels = ((Bfield, L"|B| (R,Z)"), (A_ϕ, L"A_\phi (R,Z)"),
          (A_R, L"A_R (R,Z)"), (A_Z, L"A_Z (R,Z)"))

for (n, (vals, title)) in enumerate(panels)
    ax = Axis(fig[cld(n,2), mod1(n,2)];
        xlabel = L"R", ylabel = L"Z", title = title, aspect = DataAspect())
    contour!(ax, Rgrid, Zgrid, vals; levels = 25)
    lines!(ax, boundary_R, boundary_Z; color = :red, linewidth = 3)
end

fig
```
