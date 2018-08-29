
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


function analyticA₁(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.A * sin(x[3]) + equ.C * cos(x[2])
end

function analyticA₂(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.B * sin(x[1]) + equ.A * cos(x[3])
end

function analyticA₃(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    equ.C * sin(x[2]) + equ.B * cos(x[1])
end

function analyticMetric(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    Sym[1  0  0;
        0  1  0;
        0  0  1]
end

function analyticHcoeffs(x::AbstractArray{T,1}, equ::ABC) where {T <: Number}
    Sym[1  0  0;
        0  1  0;
        0  0  1]
end
