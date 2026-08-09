
using CairoMakie
using ElectromagneticFields
using Test

# One instance of every equilibrium that has a plotting method. Loading CairoMakie
# above is what brings the Makie extension into scope, so this testset doubles as a
# check that the extension is found and precompiles.
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
    Solovev.ITER(xpoint=true),
    Solovev.NSTXdoubleX(),
)


@testset "$(rpad("Plotting",60))" begin

    for equ in plot_equilibria
        @test plot_equilibrium(equ) isa Makie.Figure
    end

    # keywords are forwarded to the axis and the contour call, and the per-panel
    # defaults must not collide with what the caller passes
    for equ in plot_equilibria
        @test plot_equilibrium(equ; levels=5, colorbar=true, colormap=:plasma) isa Makie.Figure
    end

    # drawing into an existing figure, which is how several equilibria are composed
    fig = Makie.Figure()
    for (n, equ) in enumerate(plot_equilibria)
        @test plot_equilibrium!(fig[1, n], equ) isa Makie.GridPosition
    end

    # equilibria without a plotting method
    @test_throws MethodError plot_equilibrium(PenningTrapUniform.init())

end
