
using CairoMakie
using ElectromagneticFields
using Test

# One instance of every equilibrium that has a plotting method. Loading CairoMakie above is what
# brings the extension into scope, so this testset doubles as a check that it is found at all.
const plot_equilibria = (
    ABC.init(),
    AxisymmetricTokamakCartesian.init(),
    AxisymmetricTokamakCylindrical.init(),
    AxisymmetricTokamakToroidal.init(),
    Dipole.init(),
    QuadraticPotentials.init(),
    Singular.init(),
    SolovevSymmetric.init(),
    SymmetricQuadratic.init(),
    ThetaPinch.init(),
    Solovev.ITER(),
    Solovev.NSTX(),
    Solovev.FRC(),
    Solovev.ITER(xpoint = true),
    Solovev.NSTXdoubleX()
)

# ABC is sampled on a cubic grid, so it takes neither `ny` nor `xlims`/`ylims` and has to sit out
# the tests that vary them.
const rectangular_equilibria = filter(equ -> !(equ isa ABC.ABCEquilibrium), plot_equilibria)

isaxis(x) = x isa Makie.Axis || (x isa AbstractVector && all(ax -> ax isa Makie.Axis, x))

@testset "$(rpad("Plotting",60))" begin
    for equ in plot_equilibria
        @test plot_equilibrium(equ) isa Makie.Figure
    end

    # keywords are forwarded to the axis and the contour call, and the per-panel
    # defaults must not collide with what the caller passes
    for equ in plot_equilibria
        @test plot_equilibrium(equ; levels = 5, colorbar = true, colormap = :plasma) isa
              Makie.Figure
    end

    # `size` and `figure` can both carry a figure size, so the precedence is part of the interface:
    # the size chosen for the equilibrium, then `figure`, then an explicit `size`
    figuresize(fig) = Tuple(Makie.widths(fig.scene.viewport[]))
    let equ = ThetaPinch.init()
        @test figuresize(plot_equilibrium(equ)) == (800, 400)
        @test figuresize(plot_equilibrium(equ; size = (320, 240))) == (320, 240)
        @test figuresize(plot_equilibrium(equ; figure = (; size = (360, 260)))) ==
              (360, 260)
        @test figuresize(plot_equilibrium(equ; size = (320, 240), figure = (;
            size = (360, 260)))) == (320, 240)
    end

    # A non-square grid is what guards against the value matrices being transposed:
    # Makie rejects a matrix whose dimensions do not match (length(x), length(y)),
    # while a square grid accepts a transposed one silently.
    for equ in rectangular_equilibria
        @test plot_equilibrium(equ; nx = 37, ny = 53) isa Makie.Figure
    end

    # plot ranges away from the defaults
    for equ in rectangular_equilibria
        @test plot_equilibrium(equ; xlims = (0.6, 1.4), ylims = (-0.4, +0.4)) isa
              Makie.Figure
    end

    # the singular field diverges on the z axis, so its logarithmic contour levels have
    # to cope with a grid that hits the axis, and with a window in which a component of
    # the vector potential does not change sign
    @test plot_equilibrium(Singular.init(); nx = 101, ny = 101) isa Makie.Figure
    @test plot_equilibrium(Singular.init(); nx = 101, ny = 101, colorbar = true) isa
          Makie.Figure
    @test plot_equilibrium(Singular.init(); xlims = (2.0, 3.0), ylims = (2.0, 3.0)) isa
          Makie.Figure

    # drawing into an existing figure, which is how several equilibria are composed,
    # returns the axis or axes that were created
    fig = Makie.Figure()
    for (n, equ) in enumerate(plot_equilibria)
        @test isaxis(plot_equilibrium!(fig[1, n], equ))
    end

    # every position type the extension accepts
    fig = Makie.Figure()
    @test isaxis(plot_equilibrium!(fig[1, 1], Solovev.ITER()))
    @test isaxis(plot_equilibrium!(fig[1, 2][1, 1], Solovev.NSTX()))
    @test isaxis(plot_equilibrium!(Makie.GridLayout(fig[1, 3]), ThetaPinch.init()))

    # equilibria without a plotting method
    @test_throws ArgumentError plot_equilibrium(PenningTrapUniform.init())
    @test_throws ArgumentError plot_equilibrium(PenningTrapBottle.init())
    @test_throws ArgumentError plot_equilibrium(PenningTrapAsymmetric.init())
end
