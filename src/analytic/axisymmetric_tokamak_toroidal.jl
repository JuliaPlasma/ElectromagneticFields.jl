@doc raw"""
Axisymmetric tokamak equilibrium in (r,θ,ϕ) coordinates with covariant
components of the vector potential given by
```math
A (r, \theta, \phi) = B_0 \, \bigg( 0 , \, \frac{r R_0}{\cos (\theta)} - \bigg( \frac{R_0}{\cos (\theta)} \bigg)^2 \, \ln \bigg( \frac{R}{R_0} \bigg) , \, - \frac{r^2}{2 q_0} \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (r, \theta, \phi) = \frac{B_0}{q_0} \, \bigg( 0 , \, \frac{r^2}{R}, \, q_0 R_0 \bigg)^T ,
```
where $R = R_0 + r \cos \theta$.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `q₀`: safety factor at magnetic axis
"""
module AxisymmetricTokamakToroidal

    using RecipesBase

    import ..ElectromagneticFields
    import ..ElectromagneticFields: AnalyticEquilibrium, code

    export  AxisymmetricTokamakToroidalEquilibrium

    const DEFAULT_R₀ = 1.0
    const DEFAULT_B₀ = 1.0
    const DEFAULT_q₀ = 2.0

    struct AxisymmetricTokamakToroidalEquilibrium{T <: Number} <: AnalyticEquilibrium
        name::String
        R₀::T
        B₀::T
        q₀::T

        function AxisymmetricTokamakToroidalEquilibrium{T}(R₀::T, B₀::T, q₀::T) where T <: Number
            new("AxisymmetricTokamakEquilibriumToroidal", R₀, B₀, q₀)
        end
    end

    AxisymmetricTokamakToroidalEquilibrium(R₀::T=DEFAULT_R₀, B₀::T=DEFAULT_B₀, q₀::T=DEFAULT_q₀) where T <: Number = AxisymmetricTokamakToroidalEquilibrium{T}(R₀, B₀, q₀)

    function init(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀)
        AxisymmetricTokamakToroidalEquilibrium(R₀, B₀, q₀)
    end

    function ITER()
        AxisymmetricTokamakToroidalEquilibrium(ITER_R₀, ITER_B₀, ITER_q₀)
    end

    macro code(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀)
        code(init(R₀, B₀, q₀); escape=true)
    end

    macro code_iter()
        code(ITER(); escape=true)
    end
    
    function Base.show(io::IO, equ::AxisymmetricTokamakToroidalEquilibrium)
        print(io, "Axisymmetric Tokamak Equilibrium in Toroidal Coordinates with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  q₀ = ", equ.q₀)
    end


    r(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = x[1]
    θ(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = x[2]
    ϕ(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = x[3]
    R(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = equ.R₀ + r(x,equ) * cos(θ(x,equ))
    X(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = R(x,equ) * cos(ϕ(x,equ))
    Y(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = R(x,equ) * sin(ϕ(x,equ))
    Z(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = r(x,equ) * sin(θ(x,equ))

    ElectromagneticFields.J(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = r(x,equ) * R(x,equ)
    
    ElectromagneticFields.A₁(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = + equ.B₀ * equ.R₀ * ( Z(x,equ) / R(x,equ) * cos(θ(x,equ)) - log(R(x,equ) / equ.R₀) * sin(θ(x,equ)) ) / 2
    ElectromagneticFields.A₂(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = - equ.B₀ * equ.R₀ * ( Z(x,equ) / R(x,equ) * sin(θ(x,equ)) + log(R(x,equ) / equ.R₀) * cos(θ(x,equ)) ) * r(x,equ) / 2
    ElectromagneticFields.A₃(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = + equ.B₀ * r(x,equ)^2 / equ.q₀ / 2

    ElectromagneticFields.x¹(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = X(ξ,equ)
    ElectromagneticFields.x²(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = Y(ξ,equ)
    ElectromagneticFields.x³(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = Z(ξ,equ)

    ElectromagneticFields.ξ¹(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = sqrt((sqrt(x[1]^2 + x[2]^2)-equ.R₀)^2 + x[3]^2)
    ElectromagneticFields.ξ²(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = atan(x[3], sqrt(x[1]^2 + x[2]^2)-equ.R₀)
    ElectromagneticFields.ξ³(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = atan(x[2], x[1])

    ElectromagneticFields.g₁₁(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = one(eltype(x))
    ElectromagneticFields.g₂₂(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = r(x, equ)^2
    ElectromagneticFields.g₃₃(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = R(x, equ)^2

    ElectromagneticFields.get_functions(::AxisymmetricTokamakToroidalEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ)

    function ElectromagneticFields.periodicity(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
        p = zero(x)
        p[2] = 2π
        p[3] = 2π
        return p
    end


    @recipe function f(equ::AxisymmetricTokamakToroidalEquilibrium;
                       nx = 100, ny = 120, levels = 50, size = (400,400),
                       xlims = (  0.5 * equ.R₀,   1.5 * equ.R₀),
                       ylims = (- 0.5 * equ.R₀, + 0.5 * equ.R₀))

        xgrid = LinRange(xlims[1], xlims[2], nx)
        zgrid = LinRange(ylims[1], ylims[2], ny)
        rgrid = [ξ¹(0, xgrid[i], 0.0, zgrid[j]) for i in eachindex(xgrid), j in eachindex(zgrid)]
        θgrid = [ξ²(0, xgrid[i], 0.0, zgrid[j]) for i in eachindex(xgrid), j in eachindex(zgrid)]
        pot   = [A₃(0, rgrid[i,j], θgrid[i,j], 0.0) / xgrid[i] for i in eachindex(xgrid), j in eachindex(zgrid)]

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
