"""
Axisymmetric tokamak equilibrium in (x,y,z) coordinates.

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
