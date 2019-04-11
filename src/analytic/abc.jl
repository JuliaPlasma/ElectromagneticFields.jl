@doc raw"""
ABC equilibrium in (x,y,z) coordinates with covariant components of the vector
potential given by
```math
A (x,y,z) = \big( a \, \sin(z) + c \, \cos(y) , \, b \, \sin(x) + a \, \cos(z) , \, c \, \sin(y) + b \, \cos(x) \big)^T
```
resulting in the magnetic field ``B(x,y,z) = A(x,y,z)``.

Parameters: `a`, `b`, `c`
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


function X(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    x[1]
end

function Y(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    x[2]
end

function Z(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    x[3]
end


function J(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    one(T)
end


function A₁(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.a * sin(x[3]) + equ.c * cos(x[2])
end

function A₂(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.b * sin(x[1]) + equ.a * cos(x[3])
end

function A₃(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.c * sin(x[2]) + equ.b * cos(x[1])
end


function g₁₁(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    one(T)
end

function g₂₂(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    one(T)
end

function g₃₃(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    one(T)
end


macro abc_equilibrium(a, b, c)
    generate_equilibrium_code(ABC(a, b, c); output=false)
end
