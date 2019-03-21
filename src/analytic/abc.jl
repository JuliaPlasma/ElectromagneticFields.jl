
"""
ABC equilibrium in (x,y,z) coordinates.

Parameters:
    A:
    B:
    C:
"""
struct ABC{T <: Number} <: AnalyticEquilibrium
    name::String
    A::T
    B::T
    C::T

    ABC{T}(A::T, B::T, C::T) where T <: Number = new("ABCEquilibrium", A, B, C)
end

ABC(A::T, B::T, C::T) where T <: Number = ABC{T}(A, B, C)


function Base.show(io::IO, equ::ABC)
    print(io, "ABC Equilibrium with\n")
    print(io, "  A = ", equ.A, "\n")
    print(io, "  B = ", equ.B, "\n")
    print(io, "  C = ", equ.C)
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


@inline function A₁(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.A * sin(x[3]) + equ.C * cos(x[2])
end

@inline function A₂(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.B * sin(x[1]) + equ.A * cos(x[3])
end

@inline function A₃(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.C * sin(x[2]) + equ.B * cos(x[1])
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
