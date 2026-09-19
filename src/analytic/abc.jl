
const DEFAULT_ABC_A = 1
const DEFAULT_ABC_B = 1
const DEFAULT_ABC_C = 1

@doc raw"""
Arnold-Beltrami-Childress (ABC) field in (x,y,z) coordinates with covariant components of the vector
potential given by
```math
A (x,y,z) = \big( a \, \sin(z) + c \, \cos(y) , \, b \, \sin(x) + a \, \cos(z) , \, c \, \sin(y) + b \, \cos(x) \big)^T
```
resulting in the magnetic field ``B(x,y,z) = A(x,y,z)``.

Parameters: `a`, `b`, `c`
"""
struct ABCEquilibrium{T <: Number} <: CartesianEquilibrium
    name::String
    a₀::T
    b₀::T
    c₀::T

    ABCEquilibrium{T}(a::T, b::T, c::T) where {T <: Number} = new("ABCEquilibrium", a, b, c)
end

function ABCEquilibrium(a::T = DEFAULT_ABC_A, b::T = DEFAULT_ABC_B,
        c::T = DEFAULT_ABC_C) where {T <: Number}
    ABCEquilibrium{T}(a, b, c)
end

function Base.show(io::IO, equ::ABCEquilibrium)
    print(io, "ABC Equilibrium with\n")
    print(io, "  A = ", equ.a₀, "\n")
    print(io, "  B = ", equ.b₀, "\n")
    print(io, "  C = ", equ.c₀)
end

A₁(x::AbstractVector, equ::ABCEquilibrium) = equ.a₀ * sin(x[3]) + equ.c₀ * cos(x[2])
A₂(x::AbstractVector, equ::ABCEquilibrium) = equ.b₀ * sin(x[1]) + equ.a₀ * cos(x[3])
A₃(x::AbstractVector, equ::ABCEquilibrium) = equ.c₀ * sin(x[2]) + equ.b₀ * cos(x[1])

B(x::AbstractVector, equ::ABCEquilibrium) = sqrt(A₁(x, equ)^2 + A₂(x, equ)^2 + A₃(x, equ)^2)

get_functions(::ABCEquilibrium) = (X = X, Y = Y, Z = Z)
