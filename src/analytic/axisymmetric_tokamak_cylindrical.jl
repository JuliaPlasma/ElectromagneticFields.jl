@doc raw"""
Axisymmetric tokamak equilibrium in (R,Z,ϕ) coordinates with covariant
components of the vector potential given by
```math
A (R, Z, \phi) = \frac{B_0}{2} \, \bigg( R_0 \, \frac{Z}{R} , \, - R_0 \, \ln \bigg( \frac{R}{R_0} \bigg) , \, - \frac{r^2}{q_0} \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (R, Z, \phi) = \frac{B_0}{q_0} \, \bigg( - \frac{Z}{R} , \, \frac{R - R_0}{R} , \, - q_0 R_0 \bigg)^T ,
```
where $r = \sqrt{ (R - R_0)^2 + Z^2 }$.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `q₀`: safety factor at magnetic axis
"""
module AxisymmetricTokamakCylindrical

using RecipesBase

import NaNMath: log

import ..ElectromagneticFields
import ..ElectromagneticFields: AnalyticEquilibrium, code

export AxisymmetricTokamakCylindricalEquilibrium

const DEFAULT_R₀ = 1.0
const DEFAULT_B₀ = 1.0
const DEFAULT_q₀ = 2.0

const ITER_R₀ = 6.2
const ITER_B₀ = 5.3
const ITER_q₀ = √2

struct AxisymmetricTokamakCylindricalEquilibrium{T<:Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakCylindricalEquilibrium{T}(R₀::T, B₀::T, q₀::T) where {T<:Number}
        new("AxisymmetricTokamakCylindricalEquilibrium", R₀, B₀, q₀)
    end
end

AxisymmetricTokamakCylindricalEquilibrium(R₀::T=DEFAULT_R₀, B₀::T=DEFAULT_B₀, q₀::T=DEFAULT_q₀) where {T<:Number} = AxisymmetricTokamakCylindricalEquilibrium{T}(R₀, B₀, q₀)

function init(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀)
    AxisymmetricTokamakCylindricalEquilibrium(R₀, B₀, q₀)
end

function ITER()
    AxisymmetricTokamakCylindricalEquilibrium(ITER_R₀, ITER_B₀, ITER_q₀)
end

macro code(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀)
    code(init(R₀, B₀, q₀); escape=true)
end

macro code_iter()
    code(ITER(); escape=true)
end

function Base.show(io::IO, equ::AxisymmetricTokamakCylindricalEquilibrium)
    print(io, "Axisymmetric Tokamak Equilibrium in (R,Z,ϕ) Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q₀ = ", equ.q₀)
end


R(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = x[1]
Z(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = x[2]
ϕ(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = x[3]
r²(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = (R(x, equ) - equ.R₀)^2 + Z(x, equ)^2
r(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = sqrt(r²(x, equ))
X(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = R(x, equ) * cos(ϕ(x, equ))
Y(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = R(x, equ) * sin(ϕ(x, equ))
θ(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = atan(Z(x, equ), R(x, equ) - equ.R₀)

ElectromagneticFields.J(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = R(x, equ)

ElectromagneticFields.A₁(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = +equ.B₀ * equ.R₀ * Z(x, equ) / R(x, equ) / 2
ElectromagneticFields.A₂(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = -equ.B₀ * equ.R₀ * log(R(x, equ) / equ.R₀) / 2
ElectromagneticFields.A₃(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = +equ.B₀ * r²(x, equ) / equ.q₀ / 2

ElectromagneticFields.x¹(ξ::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = X(ξ, equ)
ElectromagneticFields.x²(ξ::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = Y(ξ, equ)
ElectromagneticFields.x³(ξ::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = Z(ξ, equ)

ElectromagneticFields.ξ¹(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = sqrt(x[1]^2 + x[2]^2)
ElectromagneticFields.ξ²(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = x[3]
ElectromagneticFields.ξ³(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = atan(x[2], x[1])

ElectromagneticFields.g₁₁(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = one(eltype(x))
ElectromagneticFields.g₂₂(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = one(eltype(x))
ElectromagneticFields.g₃₃(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = R(x, equ)^2

ElectromagneticFields.get_functions(::AxisymmetricTokamakCylindricalEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)

ElectromagneticFields.minx³(ξ::AbstractVector{T}, equ::AxisymmetricTokamakCylindricalEquilibrium) where {T} = T(0)
ElectromagneticFields.maxx³(ξ::AbstractVector{T}, equ::AxisymmetricTokamakCylindricalEquilibrium) where {T} = T(2π)


@recipe function f(equ::AxisymmetricTokamakCylindricalEquilibrium;
    nx=100, ny=120, levels=50, size=(400, 400),
    xlims=(0.5 * equ.R₀, 1.5 * equ.R₀),
    ylims=(-0.5 * equ.R₀, +0.5 * equ.R₀))

    xgrid = LinRange(xlims[1], xlims[2], nx)
    zgrid = LinRange(ylims[1], ylims[2], ny)
    pot = [ElectromagneticFields.A₃([xgrid[i], zgrid[j], 0.0], equ) / xgrid[i] for i in eachindex(xgrid), j in eachindex(zgrid)]

    seriestype := :contour
    aspect_ratio := :equal
    size := size
    xlims := xlims
    ylims := ylims
    levels := levels
    legend := :none

    (xgrid, zgrid, pot')
end

end
