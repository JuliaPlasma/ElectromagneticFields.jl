module ElectromagneticFieldsMakieExt

using Makie
using Makie: AxisAspect, DataAspect, Figure, GridLayout, GridPosition, GridSubposition

using ElectromagneticFields
using ElectromagneticFields: A₁, A₂, A₃, ξ¹, ξ²
using ElectromagneticFields.ABC: ABCEquilibrium
using ElectromagneticFields.AxisymmetricTokamakCartesian: AxisymmetricTokamakCartesianEquilibrium
using ElectromagneticFields.AxisymmetricTokamakCylindrical: AxisymmetricTokamakCylindricalEquilibrium
using ElectromagneticFields.AxisymmetricTokamakToroidal: AxisymmetricTokamakToroidalEquilibrium
using ElectromagneticFields.Dipole: DipoleField
using ElectromagneticFields.QuadraticPotentials: QuadraticPotentialsField
using ElectromagneticFields.Singular: SingularEquilibrium
using ElectromagneticFields.Solovev: SolovevEquilibrium, SolovevXpointEquilibrium
using ElectromagneticFields.SolovevSymmetric: SolovevSymmetricEquilibrium
using ElectromagneticFields.SymmetricQuadratic: SymmetricQuadraticEquilibrium
using ElectromagneticFields.ThetaPinch: ThetaPinchEquilibrium

import ElectromagneticFields: plot_equilibrium, plot_equilibrium!


#
# Generic entry point
#
# Every equilibrium implements `plot_equilibrium!(position, equ; kwargs...)` and
# `figuresize(equ)`; the figure-creating method is shared by all of them.
#

const Position = Union{GridPosition,GridSubposition,GridLayout}

function plot_equilibrium(equ::ElectromagneticFields.AnalyticEquilibrium;
    size=figuresize(equ), figure=NamedTuple(), kwargs...)
    fig = Figure(; size=size, figure...)
    plot_equilibrium!(fig[1, 1], equ; kwargs...)
    fig
end


#
# Shared drawing helper
#
# Draws a single contour panel, optionally with a colorbar to its right. The
# values are expected in Makie's convention, i.e. `vals[i,j]` holds the value at
# `(xgrid[i], ygrid[j])`.
#

function contourpanel!(position::Position, xgrid, ygrid, vals;
    title="", xlabel=L"x", ylabel=L"y", levels=20,
    aspect=DataAspect(), colorbar=false, colormap=:viridis, kwargs...)

    ax = Axis(position[1, 1]; title=title, xlabel=xlabel, ylabel=ylabel, aspect=aspect)
    contour!(ax, xgrid, ygrid, vals; levels=levels, colormap=colormap, kwargs...)

    if colorbar
        # a line contour carries no single colormap that Makie could derive a
        # colorbar from, so it is built from the range of the data instead
        Colorbar(position[1, 2]; colormap=colormap, limits=extrema(filter(isfinite, vals)))
    end

    ax
end

# Merges the per-panel defaults with the keywords the caller passed, letting the
# latter win, so that e.g. `title` and `levels` can be overridden from outside
# even though every panel provides its own.
panelopts(kwargs; defaults...) = (; defaults..., kwargs...)

grid(lims, n) = LinRange(lims[1], lims[2], n)


#
# Arnold-Beltrami-Childress field
#

figuresize(::ABCEquilibrium) = (400, 1200)

function plot_equilibrium!(position::Position, equ::ABCEquilibrium;
    nx=99, ni=div(nx, 2) + 1, levels=12, kwargs...)

    lims = (0, 2π)
    xgrid = grid(lims, nx)
    Bfield = [ElectromagneticFields.ABC.B([xgrid[i], xgrid[j], xgrid[k]], equ)
              for i in eachindex(xgrid), j in eachindex(xgrid), k in eachindex(xgrid)]

    for (n, (vals, xlabel, ylabel, title)) in enumerate((
        (Bfield[:, :, ni], L"x", L"y", L"|B(x,y,\pi)|"),
        (Bfield[:, ni, :], L"x", L"z", L"|B(x,\pi,z)|"),
        (Bfield[ni, :, :], L"y", L"z", L"|B(\pi,y,z)|"),
    ))
        contourpanel!(position[n, 1], xgrid, xgrid, vals;
            panelopts(kwargs; title=title, xlabel=xlabel, ylabel=ylabel, levels=levels)...)
    end

    position
end


#
# Axisymmetric tokamak equilibria
#

figuresize(::AxisymmetricTokamakCartesianEquilibrium) = (400, 400)

function plot_equilibrium!(position::Position, equ::AxisymmetricTokamakCartesianEquilibrium;
    nx=100, ny=120, levels=50,
    xlims=(0.5 * equ.R₀, 1.5 * equ.R₀),
    ylims=(-0.5 * equ.R₀, +0.5 * equ.R₀), kwargs...)

    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₂([xgrid[i], 0.0, zgrid[j]], equ) for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel=L"x", ylabel=L"z", title=L"A_y (x,0,z)", levels=levels)...)

    position
end


figuresize(::AxisymmetricTokamakCylindricalEquilibrium) = (400, 400)

function plot_equilibrium!(position::Position, equ::AxisymmetricTokamakCylindricalEquilibrium;
    nx=100, ny=120, levels=50,
    xlims=(0.5 * equ.R₀, 1.5 * equ.R₀),
    ylims=(-0.5 * equ.R₀, +0.5 * equ.R₀), kwargs...)

    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₃([xgrid[i], zgrid[j], 0.0], equ) / xgrid[i] for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel=L"R", ylabel=L"Z", title=L"A_\phi (R,Z) / R", levels=levels)...)

    position
end


figuresize(::AxisymmetricTokamakToroidalEquilibrium) = (400, 400)

function plot_equilibrium!(position::Position, equ::AxisymmetricTokamakToroidalEquilibrium;
    nx=100, ny=120, levels=50,
    xlims=(0.5 * equ.R₀, 1.5 * equ.R₀),
    ylims=(-0.5 * equ.R₀, +0.5 * equ.R₀), kwargs...)

    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    rgrid = [ξ¹([xgrid[i], 0.0, zgrid[j]], equ) for i in eachindex(xgrid), j in eachindex(zgrid)]
    θgrid = [ξ²([xgrid[i], 0.0, zgrid[j]], equ) for i in eachindex(xgrid), j in eachindex(zgrid)]
    pot = [A₃([rgrid[i, j], θgrid[i, j], 0.0], equ) / xgrid[i] for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel=L"R", ylabel=L"Z", title=L"A_\phi (r,\theta) / R", levels=levels)...)

    position
end


#
# Dipole field
#

figuresize(::DipoleField) = (800, 400)

function plot_equilibrium!(position::Position, equ::DipoleField;
    nx=100, ny=100, levels=20,
    xlims=(-1.0, +1.0), ylims=(-1.0, +1.0), kwargs...)

    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    for (n, (component, title)) in enumerate((
        (A₁, L"A_x (x,y,1)"),
        (A₂, L"A_y (x,y,1)"),
    ))
        vals = [component([xgrid[i], ygrid[j], 1.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
        contourpanel!(position[1, n], xgrid, ygrid, vals; panelopts(kwargs; title=title, levels=levels)...)
    end

    position
end


#
# Quadratic potentials
#

figuresize(::QuadraticPotentialsField) = (1200, 400)

function plot_equilibrium!(position::Position, equ::QuadraticPotentialsField;
    nx=100, ny=100, levels=20,
    xlims=(-1.0, +1.0), ylims=(-1.0, +1.0), kwargs...)

    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    for (n, (component, title)) in enumerate((
        (A₁, L"A_x (x,y,0)"),
        (A₂, L"A_y (x,y,0)"),
        (A₃, L"A_z (x,y,0)"),
    ))
        vals = [component([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
        contourpanel!(position[1, n], xgrid, ygrid, vals; panelopts(kwargs; title=title, levels=levels)...)
    end

    position
end


#
# Singular magnetic field
#
# The vector potential and the magnetic field both diverge at the origin, so the
# contour levels are spaced logarithmically instead of linearly.
#

logrange(x1, x2, n) = collect(10^y for y in range(log10(x1), log10(x2), length=n))
doublelogrange(x1, x2, n) = vcat(-logrange(x1, x2, n), +logrange(x1, x2, n))

figuresize(::SingularEquilibrium) = (400, 1200)

function plot_equilibrium!(position::Position, equ::SingularEquilibrium;
    nx=100, ny=100, levels=25,
    xlims=(-1.0, +1.0), ylims=(-1.0, +1.0), kwargs...)

    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    pot1 = [A₁([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
    pot2 = [A₂([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
    Bfield = [ElectromagneticFields.Singular.B([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]

    for (n, (vals, title, lvls)) in enumerate((
        (pot1, L"A_x (x,y)", doublelogrange(0.1, maximum(pot1), levels)),
        (pot2, L"A_y (x,y)", doublelogrange(0.1, maximum(pot2), levels)),
        (Bfield, L"B_z (x,y)", logrange(max(0.1, minimum(Bfield)), maximum(Bfield), levels)),
    ))
        contourpanel!(position[n, 1], xgrid, ygrid, vals; panelopts(kwargs; title=title, levels=lvls)...)
    end

    position
end


#
# Symmetric quadratic magnetic field
#

figuresize(::SymmetricQuadraticEquilibrium) = (400, 1200)

function plot_equilibrium!(position::Position, equ::SymmetricQuadraticEquilibrium;
    nx=100, ny=100, levels=20,
    xlims=(-1.0, +1.0), ylims=(-1.0, +1.0), kwargs...)

    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    for (n, (component, title)) in enumerate((
        (A₁, L"A_x (x,y)"),
        (A₂, L"A_y (x,y)"),
        (ElectromagneticFields.SymmetricQuadratic.B, L"B_z (x,y)"),
    ))
        vals = [component([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
        contourpanel!(position[n, 1], xgrid, ygrid, vals; panelopts(kwargs; title=title, levels=levels)...)
    end

    position
end


#
# Theta pinch
#

figuresize(::ThetaPinchEquilibrium) = (800, 400)

function plot_equilibrium!(position::Position, equ::ThetaPinchEquilibrium;
    nx=100, ny=100, levels=20,
    xlims=(-1.0, +1.0), ylims=(-1.0, +1.0), kwargs...)

    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    for (n, (component, title)) in enumerate((
        (A₁, L"A_x (x,y)"),
        (A₂, L"A_y (x,y)"),
    ))
        vals = [component([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
        contourpanel!(position[1, n], xgrid, ygrid, vals; panelopts(kwargs; title=title, levels=levels)...)
    end

    position
end


#
# Solov'ev equilibria
#

figuresize(::SolovevEquilibrium) = (300, 400)

function plot_equilibrium!(position::Position, equ::SolovevEquilibrium;
    nx=100, ny=120, nτ=200, levels=50, boundary=true,
    xlims=(0.50, 1.50), ylims=(-0.75, +0.75), kwargs...)

    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₃([xgrid[i], zgrid[j], 0.0], equ) / xgrid[i] for i in eachindex(xgrid), j in eachindex(zgrid)]

    ax = contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel=L"R / R_0", ylabel=L"Z / R_0", levels=levels)...)

    if boundary
        # the plasma boundary is the flux surface parametrised by the inverse aspect
        # ratio ϵ, the elongation κ and the triangularity δ
        τ = LinRange(0, 2π, nτ)
        boundary_X = 1 .+ equ.ϵ .* cos.(τ .+ asin(equ.δ) .* sin.(τ))
        boundary_Y = equ.ϵ .* equ.κ .* sin.(τ)
        lines!(ax, boundary_X, boundary_Y; color=:red, linewidth=3)
    end

    position
end


figuresize(::SolovevXpointEquilibrium) = (300, 400)

function plot_equilibrium!(position::Position, equ::SolovevXpointEquilibrium;
    nx=100, ny=120, levels=50,
    xlims=(0.50, 1.50), ylims=(-0.75, +0.75), kwargs...)

    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₃([xgrid[i], zgrid[j], 0.0], equ) / xgrid[i] for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel=L"R / R_0", ylabel=L"Z / R_0", levels=levels)...)

    position
end


figuresize(::SolovevSymmetricEquilibrium) = (600, 400)

function plot_equilibrium!(position::Position, equ::SolovevSymmetricEquilibrium;
    nx=100, ny=120, levels=25,
    xlims=(equ.R₀ - 0.75, equ.R₀ + 0.75), ylims=(-0.50, +0.50), kwargs...)

    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₃([xgrid[i], zgrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; title=L"A_z (x,y)", levels=levels)...)

    position
end

end
