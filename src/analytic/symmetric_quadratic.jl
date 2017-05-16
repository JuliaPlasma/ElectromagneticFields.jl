
"""
Symmetric quadratuc equilibrium in (x,y,z) coordinates.

Parameters: none
"""
struct SymmetricQuadratic{T <: Number} <: AnalyticEquilibrium
    const name::String = "SymmetricQuadraticEquilibrium"
end


function Base.show(io::IO, equ::SymmetricQuadratic)
    print(io, "Symmetric Quadratic Equilibrium")
end


function analyticA₁(x::Vector, equ::SymmetricQuadratic)
    - 0.25 * x[2] * (2. + x[1]^2 + x[2]^2)
end

function analyticA₂(x::Vector, equ::SymmetricQuadratic)
    - 0.25 * x[1] * (2. + x[1]^2 + x[2]^2)
end

function analyticA₃(x::Vector, equ::SymmetricQuadratic)
    zero(eltype(x))
end

function analyticMetric(x::Vector, equ::SymmetricQuadratic)
    [1  0  0;
     0  1  0;
     0  0  1]
end
