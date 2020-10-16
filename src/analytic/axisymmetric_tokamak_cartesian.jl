@doc raw"""
Axisymmetric tokamak equilibrium in (x,y,z) coordinates with covariant
components of the vector potential given by
```math
A (x,y,z) = \frac{1}{2} \frac{B_0}{q_0} \, \bigg( \frac{q_0 R_0 x z - r^2 y}{R^2} , \, \frac{q_0 R_0 y z + r^2 x}{R^2} , \, - q_0 R_0 \, \ln \bigg( \frac{R}{R_0} \bigg) \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (x,y,z) = \frac{B_0}{q_0} \, \bigg( - \frac{q_0 R_0 y + x z}{R^2} , \, \frac{q_0 R_0 x - y z}{R^2} , \, \frac{R - R_0}{R} \bigg)^T ,
```
where $R = \sqrt{ x^2 + y^2 }$ and $r = \sqrt{ (R - R_0)^2 + z^2 }$.

Parameters:
* `R₀`: position of magnetic axis
* `B₀`: B-field at magnetic axis
* `q₀`: safety factor at magnetic axis
"""
module AxisymmetricTokamakCartesian

    using RecipesBase

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, code
    import ..AnalyticCartesianField: X, Y, Z

    export  AxisymmetricTokamakCartesianEquilibrium

    const DEFAULT_R₀ = 1.0
    const DEFAULT_B₀ = 1.0
    const DEFAULT_q₀ = 2.0

    const ITER_R₀ = 6.2
    const ITER_B₀ = 5.3
    const ITER_q₀ = √2

    struct AxisymmetricTokamakCartesianEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        R₀::T
        B₀::T
        q₀::T

        function AxisymmetricTokamakCartesianEquilibrium{T}(R₀::T, B₀::T, q₀::T) where T <: Number
            new("AxisymmetricTokamakCartesianEquilibrium", R₀, B₀, q₀)
        end
    end

    AxisymmetricTokamakCartesianEquilibrium(R₀::T=DEFAULT_R₀, B₀::T=DEFAULT_B₀, q₀::T=DEFAULT_q₀) where T <: Number = AxisymmetricTokamakCartesianEquilibrium{T}(R₀, B₀, q₀)

    function init(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀)
        AxisymmetricTokamakCartesianEquilibrium(R₀, B₀, q₀)
    end

    function ITER()
        AxisymmetricTokamakCartesianEquilibrium(ITER_R₀, ITER_B₀, ITER_q₀)
    end

    macro code(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀)
        code(init(R₀, B₀, q₀); escape=true)
    end

    macro code_iter()
        code(ITER(); escape=true)
    end

    function Base.show(io::IO, equ::AxisymmetricTokamakCartesianEquilibrium)
        print(io, "Axisymmetric Tokamak Equilibrium in (x,y,z) Coordinates with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  q₀ = ", equ.q₀)
    end


    R²(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = X(x,equ)^2 + Y(x,equ)^2
    r²(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = (R(x,equ) - equ.R₀)^2 + Z(x,equ)^2
    R(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = sqrt(R²(x,equ))
    r(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = sqrt(r²(x,equ))
    θ(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = atan(Z(x,equ), R(x,equ) - equ.R₀)
    ϕ(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = atan(Y(x,equ), X(x,equ))

    ElectromagneticFields.A₁(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = + equ.B₀ * (equ.R₀ * X(x,equ) * Z(x,equ) - r²(x,equ) * Y(x,equ) / equ.q₀ ) / R²(x,equ) / 2
    ElectromagneticFields.A₂(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = + equ.B₀ * (equ.R₀ * Y(x,equ) * Z(x,equ) + r²(x,equ) * X(x,equ) / equ.q₀ ) / R²(x,equ) / 2
    ElectromagneticFields.A₃(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = - equ.B₀ * equ.R₀ * log(R(x,equ) / equ.R₀) / 2

    ElectromagneticFields.get_functions(::AxisymmetricTokamakCartesianEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, R²=R², r²=r²)


    @recipe function f(equ::AxisymmetricTokamakCartesianEquilibrium;
                       nx = 100, ny = 120, levels = 50, size = (400,400),
                       xlims = (  0.5 * equ.R₀,   1.5 * equ.R₀),
                       ylims = (- 0.5 * equ.R₀, + 0.5 * equ.R₀))

        xgrid = LinRange(xlims[1], xlims[2], nx)
        zgrid = LinRange(ylims[1], ylims[2], ny)
        pot   = [ElectromagneticFields.A₂([xgrid[i], 0.0, zgrid[j]], equ) for i in eachindex(xgrid), j in eachindex(zgrid)]

        seriestype   := :contour
        aspect_ratio := :equal
        size   := size
        xlims  := xlims
        ylims  := ylims
        levels := levels
        legend := :none

        (xgrid, zgrid, pot')
    end
    
end
