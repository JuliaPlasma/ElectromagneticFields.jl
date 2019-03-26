"""
ABC equilibrium in (x,y,z) coordinates.

Parameters:
    a:
    b:
    c:
"""
struct ABC{T <: Number} <: AnalyticEquilibrium
    name::String
    a::T
    b::T
    c::T

    ABC{T}(a::T, b::T, c::T) where T <: Number = new("ABCEquilibrium", a, b, c)
end

ABC(a::T, b::T, c::T) where T <: Number = ABC{T}(a, b, c)


function Base.show(io::IO, equ::ABC)
    print(io, "ABC Equilibrium with\n")
    print(io, "  A = ", equ.a, "\n")
    print(io, "  B = ", equ.b, "\n")
    print(io, "  C = ", equ.c)
end


@inline function X(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    x[1]
end

@inline function Y(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    x[2]
end

@inline function Z(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    x[3]
end


@inline function J(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    one(T)
end


@inline function A₁(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.a * sin(x[3]) + equ.c * cos(x[2])
end

@inline function A₂(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.b * sin(x[1]) + equ.a * cos(x[3])
end

@inline function A₃(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.c * sin(x[2]) + equ.b * cos(x[1])
end


@inline function g₁₁(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    one(T)
end

@inline function g₂₂(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    one(T)
end

@inline function g₃₃(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    one(T)
end
