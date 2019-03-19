"""
θ-pinch equilibrium in (x,y,z) coordinates.

Parameters:
    B₀: B-field at magnetic axis
"""
struct ThetaPinch{T <: Number} <: AnalyticEquilibrium
    name::String
    B₀::T

    function ThetaPinch{T}(B₀::T) where T <: Number
        new("ThetaPinchEquilibrium", B₀)
    end
end

ThetaPinch(B₀::T=1.0) where T <: Number = ThetaPinch{T}(B₀)


function Base.show(io::IO, equ::ThetaPinch)
    print(io, "θ-Pinch Equilibrium in (x,y,z) Coordinates with\n")
    print(io, "  B₀ = ", equ.B₀)
end


@inline function X(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    x[1]
end

@inline function Y(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    x[2]
end

@inline function Z(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    x[3]
end

@inline function r(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    sqrt(r²(x,equ))
end

@inline function r²(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    X(x,equ)^2 + Y(x,equ)^2
end

@inline function θ(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    atan2(Y(x,equ), X(x,equ))
end


@inline function A₁(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    - equ.B₀ * Y(x,equ) / 2
end

@inline function A₂(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    + equ.B₀ * X(x,equ) / 2
end

@inline function A₃(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    zero(eltype(x))
end


@inline function g₁₁(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    one(T)
end

@inline function g₂₂(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    one(T)
end

@inline function g₃₃(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    one(T)
end


macro zpinch_equilibrium(R₀, B₀)
    generate_equilibrium_code(ThetaPinch(R₀, B₀); output=false)
end
