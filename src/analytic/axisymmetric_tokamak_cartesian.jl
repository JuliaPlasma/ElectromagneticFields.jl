
const DEFAULT_TOKAMAK_CARTESIAN_R₀ = 1
const DEFAULT_TOKAMAK_CARTESIAN_B₀ = 1
const DEFAULT_TOKAMAK_CARTESIAN_q₀ = 2

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

[`AxisymmetricTokamakCartesianITER`](@ref) returns this equilibrium with ITER's parameters.
"""
struct AxisymmetricTokamakCartesianEquilibrium{T <: Number} <: CartesianEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakCartesianEquilibrium{T}(
            R₀::T, B₀::T, q₀::T) where {T <: Number}
        new("AxisymmetricTokamakCartesianEquilibrium", R₀, B₀, q₀)
    end
end

function AxisymmetricTokamakCartesianEquilibrium(
        R₀::T = DEFAULT_TOKAMAK_CARTESIAN_R₀,
        B₀::T = DEFAULT_TOKAMAK_CARTESIAN_B₀,
        q₀::T = DEFAULT_TOKAMAK_CARTESIAN_q₀) where {T <: Number}
    AxisymmetricTokamakCartesianEquilibrium{T}(R₀, B₀, q₀)
end

"""
    AxisymmetricTokamakCartesianITER()

[`AxisymmetricTokamakCartesianEquilibrium`](@ref) with ITER's parameters, `ITER_R₀`, `ITER_B₀` and
`ITER_q₀`.
"""
function AxisymmetricTokamakCartesianITER()
    AxisymmetricTokamakCartesianEquilibrium(ITER_R₀, ITER_B₀, ITER_q₀)
end

function Base.show(io::IO, equ::AxisymmetricTokamakCartesianEquilibrium)
    print(io, "Axisymmetric Tokamak Equilibrium in (x,y,z) Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q₀ = ", equ.q₀)
end

function R²(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium)
    X(x, equ)^2 + Y(x, equ)^2
end
function r²(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium)
    (R(x, equ) - equ.R₀)^2 + Z(x, equ)^2
end
R(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = sqrt(R²(x, equ))
r(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium) = sqrt(r²(x, equ))
function θ(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium)
    atan(Z(x, equ), R(x, equ) - equ.R₀)
end
function ϕ(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium)
    atan(Y(x, equ), X(x, equ))
end

function A₁(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium)
    +equ.B₀ * (equ.R₀ * X(x, equ) * Z(x, equ) - r²(x, equ) * Y(x, equ) / equ.q₀) /
    R²(x, equ) / 2
end
function A₂(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium)
    +equ.B₀ * (equ.R₀ * Y(x, equ) * Z(x, equ) + r²(x, equ) * X(x, equ) / equ.q₀) /
    R²(x, equ) / 2
end
function A₃(x::AbstractVector, equ::AxisymmetricTokamakCartesianEquilibrium)
    # `log`, not `NaNMath.log`: unlike the cylindrical and toroidal charts this one is never
    # evaluated at R = 0, since R² = x² + y² and the chart covers all of space.
    -equ.B₀ * equ.R₀ * log(R(x, equ) / equ.R₀) / 2
end

function get_functions(::AxisymmetricTokamakCartesianEquilibrium)
    (X = X, Y = Y, Z = Z, R = R, r = r, θ = θ, ϕ = ϕ, R² = R², r² = r²)
end
