
const DEFAULT_TOKAMAK_TOROIDAL_R₀ = 1.0
const DEFAULT_TOKAMAK_TOROIDAL_B₀ = 1.0
const DEFAULT_TOKAMAK_TOROIDAL_q₀ = 2.0

@doc raw"""
Axisymmetric tokamak equilibrium in (r,θ,ϕ) coordinates with covariant
components of the vector potential given by
```math
A (r, \theta, \phi) = \frac{B_0 R_0}{2} \, \bigg( \frac{Z}{R} \cos (\theta) - \ln \bigg( \frac{R}{R_0} \bigg) \sin (\theta) , \, - r \, \bigg[ \frac{Z}{R} \sin (\theta) + \ln \bigg( \frac{R}{R_0} \bigg) \cos (\theta) \bigg] , \, + \frac{r^2}{q_0 R_0} \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (r, \theta, \phi) = \frac{B_0}{q_0} \, \bigg( 0 , \, \frac{r^2}{R}, \, q_0 R_0 \bigg)^T ,
```
where $R = R_0 + r \cos \theta$ and $Z = r \sin \theta$.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `q₀`: safety factor at magnetic axis

[`AxisymmetricTokamakToroidalITER`](@ref) returns this equilibrium with ITER's parameters.
"""
struct AxisymmetricTokamakToroidalEquilibrium{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakToroidalEquilibrium{T}(
            R₀::T, B₀::T, q₀::T) where {T <: Number}
        new("AxisymmetricTokamakEquilibriumToroidal", R₀, B₀, q₀)
    end
end

function AxisymmetricTokamakToroidalEquilibrium(
        R₀::T = DEFAULT_TOKAMAK_TOROIDAL_R₀,
        B₀::T = DEFAULT_TOKAMAK_TOROIDAL_B₀,
        q₀::T = DEFAULT_TOKAMAK_TOROIDAL_q₀) where {T <: Number}
    AxisymmetricTokamakToroidalEquilibrium{T}(R₀, B₀, q₀)
end

"""
    AxisymmetricTokamakToroidalITER()

[`AxisymmetricTokamakToroidalEquilibrium`](@ref) with ITER's parameters, `ITER_R₀`, `ITER_B₀` and
`ITER_q₀`.
"""
function AxisymmetricTokamakToroidalITER()
    AxisymmetricTokamakToroidalEquilibrium(ITER_R₀, ITER_B₀, ITER_q₀)
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
function R(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    equ.R₀ + r(x, equ) * cos(θ(x, equ))
end
function X(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    R(x, equ) * cos(ϕ(x, equ))
end
function Y(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    R(x, equ) * sin(ϕ(x, equ))
end
function Z(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    r(x, equ) * sin(θ(x, equ))
end

J(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = r(x, equ) * R(x, equ)
# (r, θ, ϕ) inherits the left-handed (R, Z, ϕ) orientation, since det ∂(R,Z)/∂(r,θ) = +r.
orientation(::AxisymmetricTokamakToroidalEquilibrium) = -1

function A₁(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    +equ.B₀ * equ.R₀ *
    (Z(x, equ) / R(x, equ) * cos(θ(x, equ)) -
     NaNMath.log(R(x, equ) / equ.R₀) * sin(θ(x, equ))) / 2
end
function A₂(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    -equ.B₀ * equ.R₀ *
    (Z(x, equ) / R(x, equ) * sin(θ(x, equ)) +
     NaNMath.log(R(x, equ) / equ.R₀) * cos(θ(x, equ))) * r(x, equ) / 2
end
function A₃(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    +equ.B₀ * r(x, equ)^2 / equ.q₀ / 2
end

x¹(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = X(ξ, equ)
x²(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = Y(ξ, equ)
x³(ξ::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = Z(ξ, equ)

function ξ¹(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    sqrt((sqrt(x[1]^2 + x[2]^2) - equ.R₀)^2 + x[3]^2)
end
function ξ²(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium)
    atan(x[3], sqrt(x[1]^2 + x[2]^2) - equ.R₀)
end
ξ³(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = atan(x[2], x[1])

g₁₁(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = one(eltype(x))
g₂₂(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = r(x, equ)^2
g₃₃(x::AbstractVector, equ::AxisymmetricTokamakToroidalEquilibrium) = R(x, equ)^2

function get_functions(::AxisymmetricTokamakToroidalEquilibrium)
    (X = X, Y = Y, Z = Z, R = R, r = r, θ = θ, ϕ = ϕ)
end

minx²(ξ::AbstractVector{T}, equ::AxisymmetricTokamakToroidalEquilibrium) where {T} = T(0)
minx³(ξ::AbstractVector{T}, equ::AxisymmetricTokamakToroidalEquilibrium) where {T} = T(0)
maxx²(ξ::AbstractVector{T}, equ::AxisymmetricTokamakToroidalEquilibrium) where {T} = T(2π)
maxx³(ξ::AbstractVector{T}, equ::AxisymmetricTokamakToroidalEquilibrium) where {T} = T(2π)

# (r, θ, ϕ): both angles are periodic. The minor radius is bounded below by zero and is not, which
# is the case that makes periodicity something other than a bounded range.
function GeometricBase.periodicity(::AxisymmetricTokamakToroidalEquilibrium)
    SVector(false, true, true)
end
