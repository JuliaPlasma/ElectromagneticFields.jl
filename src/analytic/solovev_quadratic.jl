
using SymPy: N, Sym, diff, expand, solve, subs

"""
Quadratic Solov'ev equilibrium in (R,Z,phi) coordinates.
Based on McCarthy, Physics of Plasmas 6, 3554, 1999.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    a,b: free constants
"""
struct SolovevQuadratic{T <: Number} <: AbstractSolovevEquilibrium
    name::String
    R₀::T
    B₀::T
    a::T
    b::T

    function SolovevQuadratic{T}(R₀::T, B₀::T, a::T, b::T) where T <: Number
        new("QuadraticSolovevEquilibrium", R₀, B₀, a, b)
    end
end

function SolovevQuadratic(R₀::T, B₀::T, a::T, b::T) where T <: Number
    SolovevQuadratic{T}(R₀, B₀, a, b)
end


function Base.show(io::IO, equ::SolovevQuadratic)
    print(io, "Quadratic Solovev Equilibrium with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  a  = ", equ.a, "\n")
    print(io, "  b  = ", equ.b)
end


@inline function A₁(x::AbstractArray{T,1}, equ::SolovevQuadratic) where {T <: Number}
    zero(T)
end

@inline function A₂(x::AbstractArray{T,1}, equ::SolovevQuadratic) where {T <: Number}
    zero(T)
end

@inline function A₃(x::AbstractArray{T,1}, equ::SolovevQuadratic) where {T <: Number}
    - equ.a * x[1]^4 / 8 - equ.b * x[2]^2 / 2
end


macro solovev_equilibrium_quadratic(R₀, B₀, a, b)
    generate_equilibrium_code(SolovevQuadratic(R₀, B₀, a, b); output=false)
end
