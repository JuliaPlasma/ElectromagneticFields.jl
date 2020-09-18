"""
Quadratic Solov'ev equilibrium in (R,Z,phi) coordinates.
Based on McCarthy, Physics of Plasmas 6, 3554, 1999.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `a`, b`: free constants
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


A₁(x::AbstractVector, equ::SolovevQuadratic) = zero(eltype(x))
A₂(x::AbstractVector, equ::SolovevQuadratic) = zero(eltype(x))
A₃(x::AbstractVector, equ::SolovevQuadratic) = - equ.a * R(x,equ)^4 / 8 - equ.b * Z(x,equ)^2 / 2


macro solovev_equilibrium_quadratic(R₀, B₀, a, b)
    generate_equilibrium_code(SolovevQuadratic(R₀, B₀, a, b); output=false)
end
