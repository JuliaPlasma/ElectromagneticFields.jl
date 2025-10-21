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
module AxisymmetricTokamakToroidalRegularization

import ..ElectromagneticFields
import ..ElectromagneticFields: AnalyticEquilibrium, code

export AxisymmetricTokamakToroidalRegularizationEquilibrium

const DEFAULT_R₀ = 1.0
const DEFAULT_B₀ = 1.0
const DEFAULT_q₀ = 2.0

struct AxisymmetricTokamakToroidalRegularizationEquilibrium{T<:Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakToroidalRegularizationEquilibrium{T}(R₀::T, B₀::T, q₀::T) where {T<:Number}
        new("AxisymmetricTokamakEquilibriumToroidalRegularization", R₀, B₀, q₀)
    end
end

AxisymmetricTokamakToroidalRegularizationEquilibrium(R₀::T=DEFAULT_R₀, B₀::T=DEFAULT_B₀, q₀::T=DEFAULT_q₀) where {T<:Number} = AxisymmetricTokamakToroidalRegularizationEquilibrium{T}(R₀, B₀, q₀)

function init(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀)
    AxisymmetricTokamakToroidalRegularizationEquilibrium(R₀, B₀, q₀)
end

macro code(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀)
    code(init(R₀, B₀, q₀); escape=true)
end

function Base.show(io::IO, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    print(io, "Axisymmetric Tokamak Equilibrium with Toroidal Regularization in Circular Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q₀ = ", equ.q₀)
end


r(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = x[1]
θ(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = x[2]
ϕ(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = x[3]
R(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = equ.R₀ + r(x, equ) * cos(θ(x, equ))
X(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = R(x, equ) * cos(ϕ(x, equ))
Y(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = R(x, equ) * sin(ϕ(x, equ))
Z(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = r(x, equ) * sin(θ(x, equ))

ElectromagneticFields.J(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = r(x, equ) * R(x, equ)

ElectromagneticFields.A₁(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = zero(eltype(x))
ElectromagneticFields.A₂(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = +equ.B₀ * equ.R₀ / cos(θ(x, equ))^2 * (r(x, equ) * cos(θ(x, equ)) - equ.R₀ * log(R(x, equ) / equ.R₀))
ElectromagneticFields.A₃(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = +equ.B₀ * r(x, equ)^2 / equ.q₀ / 2

ElectromagneticFields.x¹(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = X(ξ, equ)
ElectromagneticFields.x²(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = Y(ξ, equ)
ElectromagneticFields.x³(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = Z(ξ, equ)

ElectromagneticFields.ξ¹(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = sqrt((sqrt(x[1]^2 + x[2]^2) - equ.R₀)^2 + x[3]^2)
ElectromagneticFields.ξ²(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = atan(x[3], sqrt(x[1]^2 + x[2]^2) - equ.R₀)
ElectromagneticFields.ξ³(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = atan(x[2], x[1])

ElectromagneticFields.g₁₁(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = one(eltype(x))
ElectromagneticFields.g₂₂(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = r(x, equ)^2
ElectromagneticFields.g₃₃(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = R(x, equ)^2

ElectromagneticFields.get_functions(::AxisymmetricTokamakToroidalRegularizationEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ)

ElectromagneticFields.minx²(ξ::AbstractVector{T}, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) where {T} = T(0)
ElectromagneticFields.minx³(ξ::AbstractVector{T}, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) where {T} = T(0)
ElectromagneticFields.maxx²(ξ::AbstractVector{T}, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) where {T} = T(2π)
ElectromagneticFields.maxx³(ξ::AbstractVector{T}, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) where {T} = T(2π)

end
