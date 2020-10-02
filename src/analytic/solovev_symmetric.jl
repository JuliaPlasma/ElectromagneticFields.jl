"""
Symmetric Solov'ev equilibrium in (R,Z,phi) coordinates.
Based on McCarthy, Physics of Plasmas 6, 3554, 1999.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `a₀`, `b₀`: free constants
"""
module SolovevSymmetric

    using RecipesBase

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..AnalyticCartesianField: X, Y, Z

    export  SolovevSymmetricEquilibrium

    const DEFAULT_R₀ = 1.0
    const DEFAULT_B₀ = 1.0
    const DEFAULT_a₀ = 2.0
    const DEFAULT_b₀ = 0.5

    struct SolovevSymmetricEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        R₀::T
        B₀::T
        a₀::T
        b₀::T

        function SolovevSymmetricEquilibrium{T}(R₀::T, B₀::T, a₀::T, b₀::T) where T <: Number
            new("QuadraticSolovevEquilibrium", R₀, B₀, a₀, b₀)
        end
    end

    function SolovevSymmetricEquilibrium(R₀::T=DEFAULT_R₀, B₀::T=DEFAULT_B₀, a₀::T=DEFAULT_a₀, b₀::T=DEFAULT_b₀) where T <: Number
        SolovevSymmetricEquilibrium{T}(R₀, B₀, a₀, b₀)
    end


    function Base.show(io::IO, equ::SolovevSymmetricEquilibrium)
        print(io, "Quadratic Solovev Equilibrium with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  a₀ = ", equ.a₀, "\n")
        print(io, "  b₀ = ", equ.b₀)
    end

    # R²(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = X(x,equ)^2 + Y(x,equ)^2
    # r²(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = R²(x,equ)
    # R(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = sqrt(R²(x,equ))
    # r(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = sqrt(r²(x,equ))
    # θ(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = atan(Y(x,equ), X(x,equ))
    # ϕ(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = θ(x,equ)

    ElectromagneticFields.A₁(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = zero(eltype(x))
    ElectromagneticFields.A₂(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = zero(eltype(x))
    ElectromagneticFields.A₃(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = - (equ.a₀ * X(x,equ)^4 / 4 + equ.b₀ * Y(x,equ)^2 ) / 2


    macro solovev_equilibrium_quadratic(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, a₀=DEFAULT_a₀, b₀=DEFAULT_b₀)
        generate_equilibrium_code(SolovevSymmetricEquilibrium(R₀, B₀, a₀, b₀); output=false)
    end

    function init(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, a₀=DEFAULT_a₀, b₀=DEFAULT_b₀; perturbation=ZeroPerturbation())
        equilibrium = SolovevSymmetricEquilibrium(R₀, B₀, a₀, b₀)
        load_equilibrium(equilibrium, perturbation; target_module=SolovevSymmetric)
        return equilibrium
    end


    @recipe function f(equ::SolovevSymmetricEquilibrium;
                       nx = 100, ny = 120, levels = 25, size = (600,400),
                       xlims = (-0.75, +0.75),
                       ylims = (-0.50, +0.50))

        xgrid = LinRange(xlims[1], xlims[2], nx)
        zgrid = LinRange(ylims[1], ylims[2], ny)
        pot   = [A₃(0, xgrid[i], zgrid[j], 0.0) for i in eachindex(xgrid), j in eachindex(zgrid)]

        seriestype := :contour
        aspect_ratio := :equal
        size   := size
        xlims  := xlims
        ylims  := ylims
        levels := levels
        legend := :none
        xguide := "x"
        yguide := "y"

        (xgrid, zgrid, pot')
    end
    
end
