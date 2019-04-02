@doc raw"""
Symmetric quadratic equilibrium in (x,y,z) coordinates.

```math
B(x,y,z) = B_0 \, \begin{pmatrix}
0 \\
0 \\
1 + x^2 + y^2 \\
\end{pmatrix}
```

Parameters: `B₀`
"""
struct SymmetricQuadratic{T <: Number} <: AnalyticEquilibrium
    name::String
    B₀::T
    SymmetricQuadratic{T}(B₀::T) where T <: Number = new("SymmetricQuadraticEquilibrium", B₀)
end

SymmetricQuadratic(B₀::T) where T <: Number = SymmetricQuadratic{T}(B₀)


function Base.show(io::IO, equ::SymmetricQuadratic)
    print(io, "Symmetric Quadratic Equilibrium")
end


function X(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    x[1]
end

function Y(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    x[2]
end

function Z(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    x[3]
end

function R(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    r(x,equ)
end

function r(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    sqrt(r²(x,equ))
end

function r²(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    X(x,equ)^2 + Y(x,equ)^2
end

function θ(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    atan2(Y(x,equ), X(x,equ))
end

function ϕ(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    θ(x,equ)
end


function J(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    one(T)
end


function A₁(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    - equ.B₀ * x[2] * (2 + x[1]^2 + x[2]^2) / 4
end

function A₂(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    + equ.B₀ * x[1] * (2 + x[1]^2 + x[2]^2) / 4
end

function A₃(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    zero(eltype(x))
end


function φ(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
   equ.E₀ / (2π) * sin(2π * Z(x,equ))
end


function g₁₁(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    one(T)
end

function g₂₂(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    one(T)
end

function g₃₃(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    one(T)
end


macro symmetric_quadratic_equilibrium()
    generate_equilibrium_code(SymmetricQuadratic(1.); output=false)
end
