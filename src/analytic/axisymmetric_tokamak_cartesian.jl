"""
Axisymmetric tokamak equilibrium in (x,y,z) coordinates.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    q₀: safety factor at magnetic axis
"""
struct AxisymmetricTokamakCartesian{T <: Number} <: AnalyticEquilibrium
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


@inline function X(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    x[1]
end

@inline function Y(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    x[2]
end

@inline function Z(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    x[3]
end

@inline function R(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    sqrt(X(x,equ)^2 + Y(x,equ)^2)
end

@inline function r(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    sqrt(r²(x,equ))
end

@inline function r²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    (R(x,equ) - equ.R₀)^2 + Z(x,equ)^2
end

@inline function θ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    atan2(Z(x, equ), R(x, equ) - equ.R₀)
end

@inline function ϕ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    atan2(Y(x,equ), X(x,equ))
end


@inline function J(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    one(T)
end


@inline function A₁(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    + equ.B₀ * equ.R₀ * Z(x,equ) / R(x,equ) / 2
end

@inline function A₂(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    - equ.B₀ * equ.R₀ * log(R(x,equ) / equ.R₀) / 2
end

@inline function A₃(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    - equ.B₀ * r²(x,equ) / equ.q₀ / 2
end


@inline function g₁₁(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    one(T)
end

@inline function g₂₂(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    one(T)
end

@inline function g₃₃(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    one(T)
end


macro axisymmetric_tokamak_equilibrium_cartesian(R₀, B₀, q₀)
    generate_equilibrium_code(AxisymmetricTokamakCartesian(R₀, B₀, q₀); output=false)
end
