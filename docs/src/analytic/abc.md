# Arnold-Beltrami-Childress Field

```@docs
ABC
```

## Constructing the Field

The three parameters `a`, `b` and `c` all default to one. If one of them vanishes the field lines
are integrable, the flow reducing to a two-dimensional one. When all three are non-zero, regions of
chaotic field lines appear alongside regular ones — including for ``a = b = c``, the case in which
they were first studied ([Dombre et al., 1986](https://doi.org/10.1017/S0022112086002859)).

```@example abc
using CairoMakie
using ElectromagneticFields

equ = ABC.init(1.0, 0.5, 0.5)
```

## Plotting

Since the ABC field is genuinely three-dimensional, the plot shows the absolute value of the
magnetic field in the three mid-planes ``z = \pi``, ``y = \pi`` and ``x = \pi``:

```@example abc
plot_equilibrium(equ)
```

## Evaluating the Field

```@example abc
ABC.@code(1.0, 0.5, 0.5)
nothing # hide
```

Because ``B = A`` for this field, the vector potential is the more interesting quantity to look
at. Here it is sampled on a three-dimensional grid over one period:

```@example abc
nx, ny, nz = 100, 110, 120

xgrid = LinRange(0, 2π, nx)
ygrid = LinRange(0, 2π, ny)
zgrid = LinRange(0, 2π, nz)

potential(component) = [component(0.0, xgrid[i], ygrid[j], zgrid[k])
                        for i in eachindex(xgrid), j in eachindex(ygrid), k in eachindex(zgrid)]

A_x = potential(A₁)
A_y = potential(A₂)
A_z = potential(A₃)

size(A_x)
```

Cutting through the middle of the ``z`` range gives the three components in the ``(x,y)`` plane:

```@example abc
k = div(nz, 2)

fig = Figure(size = (1200, 400))

for (n, (vals, label)) in enumerate(((A_x, L"A_x"), (A_y, L"A_y"), (A_z, L"A_z")))
    ax = Axis(fig[1,n]; xlabel = L"x", ylabel = L"y", title = label, aspect = DataAspect())
    contour!(ax, xgrid, ygrid, vals[:,:,k]; levels = 12)
end

fig
```
