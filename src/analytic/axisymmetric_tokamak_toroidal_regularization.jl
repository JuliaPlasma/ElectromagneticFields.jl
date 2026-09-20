
const DEFAULT_TOKAMAK_REGULARIZATION_R₀ = 1.0
const DEFAULT_TOKAMAK_REGULARIZATION_B₀ = 1.0
const DEFAULT_TOKAMAK_REGULARIZATION_q₀ = 2.0

@doc raw"""
Axisymmetric tokamak equilibrium in (r,θ,ϕ) coordinates with covariant
components of the vector potential given by
```math
A (r, \theta, \phi) = B_0 \, \bigg( 0 , \, \bigg( \frac{R_0}{\cos (\theta)} \bigg)^2 \, \ln \bigg( \frac{R}{R_0} \bigg) - \frac{r R_0}{\cos (\theta)} , \, + \frac{r^2}{2 q_0} \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (r, \theta, \phi) = \frac{B_0}{q_0} \, \bigg( 0 , \, \frac{r^2}{R}, \, q_0 R_0 \bigg)^T ,
```
where $R = R_0 + r \cos \theta$.

This is the same magnetic field as [`AxisymmetricTokamakToroidalEquilibrium`](@ref), in a gauge
whose poloidal vector potential is regular on the magnetic axis.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `q₀`: safety factor at magnetic axis
"""
struct AxisymmetricTokamakToroidalRegularizationEquilibrium{T <: Number} <:
       AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakToroidalRegularizationEquilibrium{T}(
            R₀::T, B₀::T, q₀::T) where {T <: Number}
        new("AxisymmetricTokamakEquilibriumToroidalRegularization", R₀, B₀, q₀)
    end
end

function AxisymmetricTokamakToroidalRegularizationEquilibrium(
        R₀::T = DEFAULT_TOKAMAK_REGULARIZATION_R₀,
        B₀::T = DEFAULT_TOKAMAK_REGULARIZATION_B₀,
        q₀::T = DEFAULT_TOKAMAK_REGULARIZATION_q₀) where {T <: Number}
    AxisymmetricTokamakToroidalRegularizationEquilibrium{T}(R₀, B₀, q₀)
end

function Base.show(io::IO, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    print(io,
        "Axisymmetric Tokamak Equilibrium with Toroidal Regularization in Circular Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q₀ = ", equ.q₀)
end

r(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = x[1]
θ(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = x[2]
ϕ(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = x[3]
function R(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    equ.R₀ + r(x, equ) * cos(θ(x, equ))
end
function X(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    R(x, equ) * cos(ϕ(x, equ))
end
function Y(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    R(x, equ) * sin(ϕ(x, equ))
end
function Z(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    r(x, equ) * sin(θ(x, equ))
end

function J(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    r(x, equ) * R(x, equ)
end
# As for the unregularised toroidal chart.
orientation(::AxisymmetricTokamakToroidalRegularizationEquilibrium) = -1

function A₁(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    zero(eltype(x))
end
function A₂(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    -equ.B₀ * equ.R₀ / cos(θ(x, equ))^2 *
    (r(x, equ) * cos(θ(x, equ)) - equ.R₀ * NaNMath.log(R(x, equ) / equ.R₀))
end
function A₃(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    +equ.B₀ * r(x, equ)^2 / equ.q₀ / 2
end

x¹(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = X(ξ, equ)
x²(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = Y(ξ, equ)
x³(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) = Z(ξ, equ)

function ξ¹(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    sqrt((sqrt(x[1]^2 + x[2]^2) - equ.R₀)^2 + x[3]^2)
end
function ξ²(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    atan(x[3], sqrt(x[1]^2 + x[2]^2) - equ.R₀)
end
function ξ³(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    atan(x[2], x[1])
end

function g₁₁(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    one(eltype(x))
end
function g₂₂(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    r(x, equ)^2
end
function g₃₃(x::AbstractVector, equ::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    R(x, equ)^2
end

function get_functions(::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    (X = X, Y = Y, Z = Z, R = R, r = r, θ = θ, ϕ = ϕ)
end

function minx²(ξ::AbstractVector{T},
        equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) where {T}
    T(0)
end
function minx³(ξ::AbstractVector{T},
        equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) where {T}
    T(0)
end
function maxx²(ξ::AbstractVector{T},
        equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) where {T}
    T(2π)
end
function maxx³(ξ::AbstractVector{T},
        equ::AxisymmetricTokamakToroidalRegularizationEquilibrium) where {T}
    T(2π)
end

# (r, θ, ϕ), as the unregularised toroidal chart. The gauge changes the vector potential, not the
# coordinates.
function GeometricBase.periodicity(::AxisymmetricTokamakToroidalRegularizationEquilibrium)
    SVector(false, true, true)
end
