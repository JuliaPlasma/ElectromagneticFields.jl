module ElectromagneticFieldsMakieExt

using Makie
using Makie: AxisAspect, DataAspect, Figure, GridLayout, GridPosition, GridSubposition

using ElectromagneticFields
using ElectromagneticFields: A₁, A₂, A₃, ξ¹, ξ²
using ElectromagneticFields.ABC: ABCEquilibrium
using ElectromagneticFields.AxisymmetricTokamakCartesian:
                                                          AxisymmetricTokamakCartesianEquilibrium
using ElectromagneticFields.AxisymmetricTokamakCylindrical:
                                                            AxisymmetricTokamakCylindricalEquilibrium
using ElectromagneticFields.AxisymmetricTokamakToroidal:
                                                         AxisymmetricTokamakToroidalEquilibrium
using ElectromagneticFields.Dipole: DipoleField
using ElectromagneticFields.QuadraticPotentials: QuadraticPotentialsField
using ElectromagneticFields.Singular: SingularEquilibrium
using ElectromagneticFields.Solovev: SolovevEquilibrium, SolovevXpointEquilibrium
using ElectromagneticFields.SolovevSymmetric: SolovevSymmetricEquilibrium
using ElectromagneticFields.SymmetricQuadratic: SymmetricQuadraticEquilibrium
using ElectromagneticFields.ThetaPinch: ThetaPinchEquilibrium

import ElectromagneticFields: plot_equilibrium, plot_equilibrium!

# Every equilibrium that can be plotted implements `plot_equilibrium!(position, equ; kwargs...)` and
# `figuresize(equ)`; the figure-creating method below is shared by all of them.

const Position = Union{GridPosition, GridSubposition, GridLayout}

function plot_equilibrium(equ::ElectromagneticFields.AnalyticEquilibrium;
        size = nothing, figure = NamedTuple(), kwargs...)
    # `size` and `figure` can both carry one, so the order here is the documented precedence: the
    # size chosen for the equilibrium, then whatever `figure` holds, then an explicit `size`.
    opts = merge((; size = figuresize(equ)), figure, size === nothing ? (;) :
                                                     (; size = size))
    fig = Figure(; opts...)
    plot_equilibrium!(fig[1, 1], equ; kwargs...)
    fig
end

# Fields without a plotting method fail on `figuresize` before they ever reach
# `plot_equilibrium!`, so that is where the explanation belongs.
function figuresize(equ::ElectromagneticFields.AnalyticEquilibrium)
    throw(ArgumentError(
        "no plotting method is defined for $(typeof(equ)); see the Plotting page of the " *
        "documentation for how to sample and plot such a field directly"))
end

# Shared drawing helper. The values are expected in Makie's convention, i.e. `vals[i,j]` holds the
# value at `(xgrid[i], ygrid[j])`.

function contourpanel!(position::Position, xgrid, ygrid, vals;
        title = "", xlabel = L"x", ylabel = L"y", levels = 20,
        aspect = DataAspect(), colorbar = false, colormap = :viridis, kwargs...)
    ax = Axis(
        position[1, 1]; title = title, xlabel = xlabel, ylabel = ylabel, aspect = aspect)
    contour!(ax, xgrid, ygrid, vals; levels = levels, colormap = colormap, kwargs...)

    if colorbar
        # a line contour carries no colormap Makie could derive a colorbar from, so the range of
        # the data stands in for it
        finite = filter(isfinite, vec(vals))
        if !isempty(finite)
            lo, hi = extrema(finite)
            # Makie rejects coinciding limits, which a panel constant over the grid produces
            lo == hi && ((lo, hi) = (lo - one(lo), hi + one(hi)))
            Colorbar(position[1, 2]; colormap = colormap, limits = (lo, hi))
        end
    end

    ax
end

# Caller keywords win over the per-panel defaults, so e.g. `title` and `levels` can be overridden
# even though every panel provides its own.
panelopts(kwargs; defaults...) = (; defaults..., kwargs...)

grid(lims, n) = LinRange(lims[1], lims[2], n)

# Arnold-Beltrami-Childress field

figuresize(::ABCEquilibrium) = (400, 1200)

function plot_equilibrium!(position::Position, equ::ABCEquilibrium;
        nx = 99, ni = div(nx, 2) + 1, levels = 12, kwargs...)
    lims = (0, 2π)
    xgrid = grid(lims, nx)

    # Only the three mid-planes are shown, so only their 3nx² values are evaluated rather than the
    # full nx³ cube. `ni` picks the grid point closest to π, which for odd `nx` is π exactly.
    B(x, y, z) = ElectromagneticFields.ABC.B([x, y, z], equ)

    Bxy = [B(xgrid[i], xgrid[j], xgrid[ni])
           for i in eachindex(xgrid), j in eachindex(xgrid)]
    Bxz = [B(xgrid[i], xgrid[ni], xgrid[k])
           for i in eachindex(xgrid), k in eachindex(xgrid)]
    Byz = [B(xgrid[ni], xgrid[j], xgrid[k])
           for j in eachindex(xgrid), k in eachindex(xgrid)]

    axs = Axis[]

    for (n, (vals, xlabel, ylabel, title)) in enumerate((
        (Bxy, L"x", L"y", L"|B(x,y,\pi)|"),
        (Bxz, L"x", L"z", L"|B(x,\pi,z)|"),
        (Byz, L"y", L"z", L"|B(\pi,y,z)|")
    ))
        push!(axs,
            contourpanel!(position[n, 1], xgrid, xgrid, vals;
                panelopts(kwargs; title = title, xlabel = xlabel,
                    ylabel = ylabel, levels = levels)...))
    end

    axs
end

# Axisymmetric tokamak equilibria

figuresize(::AxisymmetricTokamakCartesianEquilibrium) = (400, 400)

function plot_equilibrium!(
        position::Position, equ::AxisymmetricTokamakCartesianEquilibrium;
        nx = 100, ny = 120, levels = 50,
        xlims = (0.5 * equ.R₀, 1.5 * equ.R₀),
        ylims = (-0.5 * equ.R₀, +0.5 * equ.R₀), kwargs...)
    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    # at y = 0 the cartesian A_y is the physical toroidal component, so multiplying by R gives
    # the covariant one, which is the flux function ψ whose contours are the flux surfaces
    pot = [A₂([xgrid[i], 0.0, zgrid[j]], equ) * xgrid[i]
           for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel = L"x", ylabel = L"z",
            title = L"R \, A_y (x,0,z)", levels = levels)...)
end

figuresize(::AxisymmetricTokamakCylindricalEquilibrium) = (400, 400)

function plot_equilibrium!(
        position::Position, equ::AxisymmetricTokamakCylindricalEquilibrium;
        nx = 100, ny = 120, levels = 50,
        xlims = (0.5 * equ.R₀, 1.5 * equ.R₀),
        ylims = (-0.5 * equ.R₀, +0.5 * equ.R₀), kwargs...)
    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₃([xgrid[i], zgrid[j], 0.0], equ)
           for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel = L"R", ylabel = L"Z",
            title = L"A_\phi (R,Z)", levels = levels)...)
end

figuresize(::AxisymmetricTokamakToroidalEquilibrium) = (400, 400)

function plot_equilibrium!(position::Position, equ::AxisymmetricTokamakToroidalEquilibrium;
        nx = 100, ny = 120, levels = 50,
        xlims = (0.5 * equ.R₀, 1.5 * equ.R₀),
        ylims = (-0.5 * equ.R₀, +0.5 * equ.R₀), kwargs...)
    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    rgrid = [ξ¹([xgrid[i], 0.0, zgrid[j]], equ)
             for i in eachindex(xgrid), j in eachindex(zgrid)]
    θgrid = [ξ²([xgrid[i], 0.0, zgrid[j]], equ)
             for i in eachindex(xgrid), j in eachindex(zgrid)]
    pot = [A₃([rgrid[i, j], θgrid[i, j], 0.0], equ)
           for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel = L"R", ylabel = L"Z",
            title = L"A_\phi (r,\theta)", levels = levels)...)
end

# Dipole field

figuresize(::DipoleField) = (800, 400)

function plot_equilibrium!(position::Position, equ::DipoleField;
        nx = 100, ny = 100, levels = 20,
        xlims = (-1.0, +1.0), ylims = (-1.0, +1.0), kwargs...)
    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    axs = Axis[]

    for (n, (component, title)) in enumerate((
        (A₁, L"A_x (x,y,1)"),
        (A₂, L"A_y (x,y,1)")
    ))
        vals = [component([xgrid[i], ygrid[j], 1.0], equ)
                for i in eachindex(xgrid), j in eachindex(ygrid)]
        push!(axs,
            contourpanel!(position[1, n], xgrid, ygrid, vals;
                panelopts(kwargs; title = title, levels = levels)...))
    end

    axs
end

# Quadratic potentials

figuresize(::QuadraticPotentialsField) = (1200, 400)

function plot_equilibrium!(position::Position, equ::QuadraticPotentialsField;
        nx = 100, ny = 100, levels = 20,
        xlims = (-1.0, +1.0), ylims = (-1.0, +1.0), kwargs...)
    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    axs = Axis[]

    for (n, (component, title)) in enumerate((
        (A₁, L"A_x (x,y,0)"),
        (A₂, L"A_y (x,y,0)"),
        (A₃, L"A_z (x,y,0)")
    ))
        vals = [component([xgrid[i], ygrid[j], 0.0], equ)
                for i in eachindex(xgrid), j in eachindex(ygrid)]
        push!(axs,
            contourpanel!(position[1, n], xgrid, ygrid, vals;
                panelopts(kwargs; title = title, levels = levels)...))
    end

    axs
end

# Singular magnetic field
#
# The vector potential and the magnetic field both diverge at the origin, so the contour levels are
# spaced logarithmically instead of linearly.

logrange(x1, x2, n) = collect(10^y for y in range(log10(x1), log10(x2), length = n))
doublelogrange(x1, x2, n) = sort!(vcat(-logrange(x1, x2, n), +logrange(x1, x2, n)))

# Taken over the finite magnitudes rather than over `maximum` directly, which copes both with a grid
# that meets the singular line, where the values are Inf or NaN, and with a component that keeps one
# sign over the plot range and so has a negative maximum. The lower bound keeps the range
# non-degenerate.
function logextent(vals; lo = 0.1)
    hi = maximum(abs, Iterators.filter(isfinite, vals); init = lo)
    max(hi, 10 * lo)
end

figuresize(::SingularEquilibrium) = (400, 1200)

function plot_equilibrium!(position::Position, equ::SingularEquilibrium;
        nx = 100, ny = 100, levels = 25,
        xlims = (-1.0, +1.0), ylims = (-1.0, +1.0), kwargs...)
    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    pot1 = [A₁([xgrid[i], ygrid[j], 0.0], equ)
            for i in eachindex(xgrid), j in eachindex(ygrid)]
    pot2 = [A₂([xgrid[i], ygrid[j], 0.0], equ)
            for i in eachindex(xgrid), j in eachindex(ygrid)]
    Bfield = [ElectromagneticFields.Singular.B([xgrid[i], ygrid[j], 0.0], equ)
              for i in eachindex(xgrid), j in eachindex(ygrid)]

    Blo = max(0.1, minimum(abs, Iterators.filter(isfinite, Bfield); init = 0.1))

    axs = Axis[]

    for (n, (vals, title, lvls)) in enumerate((
        (pot1, L"A_x (x,y)", doublelogrange(0.1, logextent(pot1), levels)),
        (pot2, L"A_y (x,y)", doublelogrange(0.1, logextent(pot2), levels)),
        (Bfield, L"B_z (x,y)", logrange(Blo, logextent(Bfield; lo = Blo), levels))
    ))
        push!(axs,
            contourpanel!(position[n, 1], xgrid, ygrid, vals;
                panelopts(kwargs; title = title, levels = lvls)...))
    end

    axs
end

# Symmetric quadratic magnetic field

figuresize(::SymmetricQuadraticEquilibrium) = (400, 1200)

function plot_equilibrium!(position::Position, equ::SymmetricQuadraticEquilibrium;
        nx = 100, ny = 100, levels = 20,
        xlims = (-1.0, +1.0), ylims = (-1.0, +1.0), kwargs...)
    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    axs = Axis[]

    for (n, (component, title)) in enumerate((
        (A₁, L"A_x (x,y)"),
        (A₂, L"A_y (x,y)"),
        (ElectromagneticFields.SymmetricQuadratic.B, L"B_z (x,y)")
    ))
        vals = [component([xgrid[i], ygrid[j], 0.0], equ)
                for i in eachindex(xgrid), j in eachindex(ygrid)]
        push!(axs,
            contourpanel!(position[n, 1], xgrid, ygrid, vals;
                panelopts(kwargs; title = title, levels = levels)...))
    end

    axs
end

# Theta pinch

figuresize(::ThetaPinchEquilibrium) = (800, 400)

function plot_equilibrium!(position::Position, equ::ThetaPinchEquilibrium;
        nx = 100, ny = 100, levels = 20,
        xlims = (-1.0, +1.0), ylims = (-1.0, +1.0), kwargs...)
    xgrid = grid(xlims, nx)
    ygrid = grid(ylims, ny)

    axs = Axis[]

    for (n, (component, title)) in enumerate((
        (A₁, L"A_x (x,y)"),
        (A₂, L"A_y (x,y)")
    ))
        vals = [component([xgrid[i], ygrid[j], 0.0], equ)
                for i in eachindex(xgrid), j in eachindex(ygrid)]
        push!(axs,
            contourpanel!(position[1, n], xgrid, ygrid, vals;
                panelopts(kwargs; title = title, levels = levels)...))
    end

    axs
end

# Solov'ev equilibria
#
# `ψ` vanishes on the plasma boundary and is negative inside it, reaching its minimum on the
# magnetic axis, while outside it grows without any bound the plot window imposes. Spreading `n`
# levels evenly over the sampled range therefore spends nearly all of them on the far field and
# leaves the flux surfaces of the plasma to a handful, so the spacing is anchored to the flux at the
# axis instead: uniform in `ψ`, one level exactly on the boundary, and `n / (1 + OUTER_FLUX_SPAN)`
# of them inside it. An explicit vector of levels overrides this, as it does everywhere else.

const OUTER_FLUX_SPAN = 3

function fluxlevels(vals, n::Integer)
    lo, hi = extrema(Iterators.filter(isfinite, vals))
    lo < 0 < hi || return n
    nin = max(1, round(Int, (n - 1) / (1 + OUTER_FLUX_SPAN)))
    Δψ = -lo / nin
    Δψ .* ((-nin):(n - 1 - nin))
end

fluxlevels(_, levels) = levels

figuresize(::SolovevEquilibrium) = (300, 400)

function plot_equilibrium!(position::Position, equ::SolovevEquilibrium;
        nx = 100, ny = 120, nτ = 200, levels = 40, boundary = true,
        xlims = (0.50, 1.50), ylims = (-0.75, +0.75), kwargs...)
    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₃([xgrid[i], zgrid[j], 0.0], equ)
           for i in eachindex(xgrid), j in eachindex(zgrid)]

    ax = contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel = L"R / R_0", ylabel = L"Z / R_0", levels = fluxlevels(pot, levels))...)

    if boundary
        # the plasma boundary is the flux surface parametrised by the inverse aspect
        # ratio ϵ, the elongation κ and the triangularity δ
        τ = LinRange(0, 2π, nτ)
        boundary_X = 1 .+ equ.ϵ .* cos.(τ .+ asin(equ.δ) .* sin.(τ))
        boundary_Y = equ.ϵ .* equ.κ .* sin.(τ)
        lines!(ax, boundary_X, boundary_Y; color = :red, linewidth = 3)
    end

    ax
end

figuresize(::SolovevXpointEquilibrium) = (300, 400)

function plot_equilibrium!(position::Position, equ::SolovevXpointEquilibrium;
        nx = 100, ny = 120, levels = 40,
        xlims = (0.50, 1.50), ylims = (-0.75, +0.75), kwargs...)
    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₃([xgrid[i], zgrid[j], 0.0], equ)
           for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; xlabel = L"R / R_0", ylabel = L"Z / R_0", levels = fluxlevels(pot, levels))...)
end

figuresize(::SolovevSymmetricEquilibrium) = (600, 400)

function plot_equilibrium!(position::Position, equ::SolovevSymmetricEquilibrium;
        nx = 100, ny = 120, levels = 25,
        # the flux function is quartic in `R₀ + x`, so the magnetic axis sits at `x = -R₀`
        xlims = (-equ.R₀ - 0.75, -equ.R₀ + 0.75), ylims = (-0.50, +0.50), kwargs...)
    xgrid = grid(xlims, nx)
    zgrid = grid(ylims, ny)
    pot = [A₃([xgrid[i], zgrid[j], 0.0], equ)
           for i in eachindex(xgrid), j in eachindex(zgrid)]

    contourpanel!(position[1, 1], xgrid, zgrid, pot;
        panelopts(kwargs; title = L"A_z (x,y)", levels = levels)...)
end

end
