
"""
ABC equilibrium in (x,y,z) coordinates.

Parameters:
    A:
    B:
    C:
"""
struct ABC{T <: Number} <: AnalyticEquilibrium
    const name::String = "ABCEquilibrium"
    A::T
    B::T
    C::T

    ABC{T}(A::T, B::T, C::T) where T <: Number = new(A, B, C)
end

ABC(A::T, B::T, C::T) where T <: Number = ABC{T}(A, B, C)


function Base.show(io::IO, equ::ABC)
    print(io, "ABC Equilibrium with\n")
    print(io, "  A = ", equ.A, "\n")
    print(io, "  B = ", equ.B, "\n")
    print(io, "  C = ", equ.C)
end


function analyticA₁(x::Vector, equ::ABC)
    equ.A * sin(x[3]) + equ.C * cos(x[2])
end

function analyticA₂(x::Vector, equ::ABC)
    equ.B * sin(x[1]) + equ.A * cos(x[3])
end

function analyticA₃(x::Vector, equ::ABC)
    equ.C * sin(x[2]) + equ.B * cos(x[1])
end

function analyticMetric(x::Vector, equ::ABC)
    [1  0  0;
     0  1  0;
     0  0  1]
end
