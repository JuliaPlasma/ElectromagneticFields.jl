
using SymPy: Sym

"""
Symmetric quadratic equilibrium in (x,y,z) coordinates.

```math
B(x,y,z) = B_0 \\, \\begin{pmatrix}
0 \\\\
0 \\\\
1 + x^2 + y^2 \\\\
\\end{pmatrix}
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


function analyticA₁(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    - equ.B₀ * x[2] * (2 + x[1]^2 + x[2]^2) / 4
end

function analyticA₂(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    + equ.B₀ * x[1] * (2 + x[1]^2 + x[2]^2) / 4
end

function analyticA₃(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    zero(eltype(x))
end

function analyticMetric(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    Sym[1  0  0;
        0  1  0;
        0  0  1]
end

function analyticHcoeffs(x::AbstractArray{T,1}, equ::SymmetricQuadratic) where {T <: Number}
    Sym[1  0  0;
        0  1  0;
        0  0  1]
end

macro symmetric_quadratic_equilibrium()
    generate_equilibrium_code(SymmetricQuadratic(1.); output=false)
end
