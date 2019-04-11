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
struct AxisymmetricTokamakCartesian{T <: Number} <: CartesianEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakCartesian{T}(R₀::T, B₀::T, q₀::T) where T <: Number
        new("AxisymmetricTokamakCartesianEquilibrium", R₀, B₀, q₀)
    end
end

AxisymmetricTokamakCartesian(R₀::T=1.0, B₀::T=1.0, q₀::T=2.0) where T <: Number = AxisymmetricTokamakCartesian{T}(R₀, B₀, q₀)


function Base.show(io::IO, equ::AxisymmetricTokamakCartesian)
    print(io, "Axisymmetric Tokamak Equilibrium in (x,y,z) Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q₀ = ", equ.q₀)
end


function R(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    sqrt(R²(x,equ))
end

function R²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    X(x,equ)^2 + Y(x,equ)^2
end

function r(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    sqrt(r²(x,equ))
end

function r²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    (R(x,equ) - equ.R₀)^2 + Z(x,equ)^2
end

function θ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    atan2(Z(x, equ), R(x, equ) - equ.R₀)
end

function ϕ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    atan2(Y(x,equ), X(x,equ))
end


function A₁(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    + equ.B₀ * (equ.R₀ * X(x,equ) * Z(x,equ) - r²(x,equ) * Y(x,equ) / equ.q₀ ) / R²(x,equ) / 2
end

function A₂(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    + equ.B₀ * (equ.R₀ * Y(x,equ) * Z(x,equ) + r²(x,equ) * X(x,equ) / equ.q₀ ) / R²(x,equ) / 2
end

function A₃(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    - equ.B₀ * equ.R₀ * log(R(x,equ) / equ.R₀) / 2
end


macro axisymmetric_tokamak_equilibrium_cartesian(R₀, B₀, q₀)
    generate_equilibrium_code(AxisymmetricTokamakCartesian(R₀, B₀, q₀); output=false)
end
