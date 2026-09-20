
const DEFAULT_TOKAMAK_CYLINDRICAL_R₀ = 1.0
const DEFAULT_TOKAMAK_CYLINDRICAL_B₀ = 1.0
const DEFAULT_TOKAMAK_CYLINDRICAL_q₀ = 2.0

@doc raw"""
Axisymmetric tokamak equilibrium in (R,Z,ϕ) coordinates with covariant
components of the vector potential given by
```math
A (R, Z, \phi) = \frac{B_0}{2} \, \bigg( R_0 \, \frac{Z}{R} , \, - R_0 \, \ln \bigg( \frac{R}{R_0} \bigg) , \, + \frac{r^2}{q_0} \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (R, Z, \phi) = \frac{B_0}{q_0} \, \bigg( - \frac{Z}{R} , \, \frac{R - R_0}{R} , \, + q_0 R_0 \bigg)^T ,
```
where $r = \sqrt{ (R - R_0)^2 + Z^2 }$.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `q₀`: safety factor at magnetic axis

[`AxisymmetricTokamakCylindricalITER`](@ref) returns this equilibrium with ITER's parameters.
"""
struct AxisymmetricTokamakCylindricalEquilibrium{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakCylindricalEquilibrium{T}(
            R₀::T, B₀::T, q₀::T) where {T <: Number}
        new("AxisymmetricTokamakCylindricalEquilibrium", R₀, B₀, q₀)
    end
end

function AxisymmetricTokamakCylindricalEquilibrium(
        R₀::T = DEFAULT_TOKAMAK_CYLINDRICAL_R₀,
        B₀::T = DEFAULT_TOKAMAK_CYLINDRICAL_B₀,
        q₀::T = DEFAULT_TOKAMAK_CYLINDRICAL_q₀) where {T <: Number}
    AxisymmetricTokamakCylindricalEquilibrium{T}(R₀, B₀, q₀)
end

"""
    AxisymmetricTokamakCylindricalITER()

[`AxisymmetricTokamakCylindricalEquilibrium`](@ref) with ITER's parameters, `ITER_R₀`, `ITER_B₀`
and `ITER_q₀`.
"""
function AxisymmetricTokamakCylindricalITER()
    AxisymmetricTokamakCylindricalEquilibrium(ITER_R₀, ITER_B₀, ITER_q₀)
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
function r²(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium)
    (R(x, equ) - equ.R₀)^2 + Z(x, equ)^2
end
r(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = sqrt(r²(x, equ))
function X(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium)
    R(x, equ) * cos(ϕ(x, equ))
end
function Y(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium)
    R(x, equ) * sin(ϕ(x, equ))
end
function θ(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium)
    atan(Z(x, equ), R(x, equ) - equ.R₀)
end

J(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = R(x, equ)
# (R, Z, ϕ) is left-handed; the right-handed ordering would be (R, ϕ, Z). See `orientation`.
orientation(::AxisymmetricTokamakCylindricalEquilibrium) = -1

function A₁(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium)
    +equ.B₀ * equ.R₀ * Z(x, equ) / R(x, equ) / 2
end
function A₂(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium)
    -equ.B₀ * equ.R₀ * NaNMath.log(R(x, equ) / equ.R₀) / 2
end
function A₃(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium)
    +equ.B₀ * r²(x, equ) / equ.q₀ / 2
end

x¹(ξ::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = X(ξ, equ)
x²(ξ::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = Y(ξ, equ)
x³(ξ::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = Z(ξ, equ)

function ξ¹(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium)
    sqrt(x[1]^2 + x[2]^2)
end
ξ²(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = x[3]
ξ³(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = atan(x[2], x[1])

g₁₁(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = one(eltype(x))
g₂₂(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = one(eltype(x))
g₃₃(x::AbstractVector, equ::AxisymmetricTokamakCylindricalEquilibrium) = R(x, equ)^2

function get_functions(::AxisymmetricTokamakCylindricalEquilibrium)
    (X = X, Y = Y, Z = Z, R = R, r = r, θ = θ, ϕ = ϕ, r² = r²)
end

function minx³(ξ::AbstractVector{T},
        equ::AxisymmetricTokamakCylindricalEquilibrium) where {T}
    T(0)
end
function maxx³(ξ::AbstractVector{T},
        equ::AxisymmetricTokamakCylindricalEquilibrium) where {T}
    T(2π)
end

# (R, Z, ϕ): the toroidal angle is periodic, the two poloidal coordinates are not.
function GeometricBase.periodicity(::AxisymmetricTokamakCylindricalEquilibrium)
    SVector(false, false, true)
end
