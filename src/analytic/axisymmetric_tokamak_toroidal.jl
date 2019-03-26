"""
Axisymmetric tokamak equilibrium in (r,θ,ϕ) coordinates.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    q₀: safety factor at magnetic axis
"""
struct AxisymmetricTokamakToroidal{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakToroidal{T}(R₀::T, B₀::T, q₀::T) where T <: Number
        new("AxisymmetricTokamakEquilibriumToroidal", R₀, B₀, q₀)
    end
end

AxisymmetricTokamakToroidal(R₀::T=1.0, B₀::T=1.0, q₀::T=2.0) where T <: Number = AxisymmetricTokamakToroidal{T}(R₀, B₀, q₀)


function Base.show(io::IO, equ::AxisymmetricTokamakToroidal)
    print(io, "Axisymmetric Tokamak Equilibrium in Toroidal Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q₀ = ", equ.q₀)
end


@inline function X(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    R(x,equ) * cos(ϕ(x,equ))
end

@inline function Y(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    R(x,equ) * sin(ϕ(x,equ))
end

@inline function Z(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    r(x,equ) * sin(θ(x,equ))
end

@inline function R(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    equ.R₀ + r(x,equ) * cos(θ(x,equ))
end

@inline function R²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    R(x,equ)^2
end

@inline function r(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    x[1]
end

@inline function r²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    r(x,equ)^2
end

@inline function θ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    x[2]
end

@inline function ϕ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    x[3]
end


@inline function J(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    r(x,equ) * R(x,equ)
end


@inline function A₁(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    zero(T)
end

@inline function A₂(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    equ.B₀ * equ.R₀ / cos(θ(x,equ))^2 * ( r(x,equ) * cos(θ(x,equ)) - equ.R₀ * log(R(x,equ) / equ.R₀) )
end

@inline function A₃(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    - equ.B₀ * r²(x,equ) / equ.q₀ / 2
end


@inline function g₁₁(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    one(T)
end

@inline function g₂₂(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    r²(x, equ)
end

@inline function g₃₃(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    R²(x, equ)
end


macro axisymmetric_tokamak_equilibrium_toroidal(R₀, B₀, q₀)
    generate_equilibrium_code(AxisymmetricTokamakToroidal(R₀, B₀, q₀); output=false)
end
