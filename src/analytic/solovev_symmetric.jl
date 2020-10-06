@doc raw"""
Symmetric Solov'ev equilibrium in (R,Z,phi) coordinates.
Based on McCarthy, Physics of Plasmas 6, 3554, 1999.

The covariant components of the vector potential are given by
```math
A (x, y) = \frac{B_0}{2} \, \bigg( 0 , \, 0 , \, - \frac{\alpha}{4} (R_0 + x)^4 - \beta y^2 \bigg)^T ,
```

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `α`, `β`: free constants
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
    const DEFAULT_α  = 2.0
    const DEFAULT_β  = 0.5

    struct SolovevSymmetricEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        R₀::T
        B₀::T
        α::T
        β::T

        function SolovevSymmetricEquilibrium{T}(R₀::T, B₀::T, α::T, β::T) where T <: Number
            new("QuadraticSolovevEquilibrium", R₀, B₀, α, β)
        end
    end

    function SolovevSymmetricEquilibrium(R₀::T=DEFAULT_R₀, B₀::T=DEFAULT_B₀, α::T=DEFAULT_α, β::T=DEFAULT_β) where T <: Number
        SolovevSymmetricEquilibrium{T}(R₀, B₀, α, β)
    end


    function Base.show(io::IO, equ::SolovevSymmetricEquilibrium)
        print(io, "Quadratic Solovev Equilibrium with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  α  = ", equ.α,  "\n")
        print(io, "  β  = ", equ.β)
    end

    # R²(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = X(x,equ)^2 + Y(x,equ)^2
    # r²(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = R²(x,equ)
    # R(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = sqrt(R²(x,equ))
    # r(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = sqrt(r²(x,equ))
    # θ(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = atan(Y(x,equ), X(x,equ))
    # ϕ(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = θ(x,equ)

    ElectromagneticFields.A₁(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = zero(eltype(x))
    ElectromagneticFields.A₂(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = zero(eltype(x))
    ElectromagneticFields.A₃(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = - equ.B₀ * (equ.α * (equ.R₀ + X(x,equ))^4 / 4 + equ.β * Y(x,equ)^2 ) / 2


    macro solovev_equilibrium_quadratic(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, α=DEFAULT_α, β=DEFAULT_β)
        generate_equilibrium_code(SolovevSymmetricEquilibrium(R₀, B₀, α, β); output=false)
    end

    function init(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, α=DEFAULT_α, β=DEFAULT_β; perturbation=ZeroPerturbation())
        equilibrium = SolovevSymmetricEquilibrium(R₀, B₀, α, β)
        load_equilibrium(equilibrium, perturbation; target_module=SolovevSymmetric)
        return equilibrium
    end


    @recipe function f(equ::SolovevSymmetricEquilibrium;
                       nx = 100, ny = 120, levels = 25, size = (600,400),
                       xlims = (equ.R₀-0.75, equ.R₀+0.75),
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
