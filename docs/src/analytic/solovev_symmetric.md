# Symmetric Solov'ev Equilibrium

```@docs
SolovevSymmetric
```

## Constructing the Field

The parameters are the position of the magnetic axis ``R_0``, the field strength ``B_0``, and the
two coefficients ``\alpha`` and ``\beta`` of the flux function:

```@example solsym
using CairoMakie
using ElectromagneticFields

equ = SolovevSymmetric.init(0.0, 1.0, 2.0, 0.5)
```

## Plotting

Unlike the [general Solov'ev equilibrium](@ref "Solov'ev Equilibrium"), this one is written in
cartesian coordinates and is symmetric both up/down and left/right. Its flux surfaces are the
contours of ``A_z``:

```@example solsym
plot_equilibrium(equ)
```

``\alpha`` weights the quartic dependence on ``x``, ``\beta`` the quadratic one on ``y``, so it is
their ratio that decides the shape of the surfaces — the larger ``\beta / \alpha``, the flatter
they get:

```@example solsym
fig = Figure(size = (1200, 400))

for (n, (α, β)) in enumerate(((2.0, 0.5), (1.0, 1.0), (0.5, 2.0)))
    plot_equilibrium!(fig[1,n], SolovevSymmetric.init(0.0, 1.0, α, β);
        title = "α = $α, β = $β")
end

fig
```
