"""
θ-pinch equilibrium in (x,y,z) coordinates.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
"""
struct ThetaPinch{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T

    function ThetaPinch{T}(R₀::T, B₀::T) where T <: Number
        new("ThetaPinchEquilibrium", R₀, B₀)
    end
end

ThetaPinch(R₀::T=0.0, B₀::T=1.0) where T <: Number = ThetaPinch{T}(R₀, B₀)


function Base.show(io::IO, equ::ThetaPinch)
    print(io, "θ-Pinch Equilibrium in (x,y,z) Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀)
end


function r²(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    (x[1] - equ.R₀)^2 + x[2]^2
end

function r(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    sqrt(r²(x, equ))
end

function θ(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    atan2(x[2], x[1] - equ.R₀)
end

function R(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    x[1]
end

function Z(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    x[2]
end

function ϕ(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    x[3]
    # TODO ! wrong !
end


function analyticA₁(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    - equ.B₀ * x[2] / 2
end

function analyticA₂(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    + equ.B₀ * x[1] / 2
end

function analyticA₃(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    zero(eltype(x))
end

function analyticMetric(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    Sym[1  0  0;
        0  1  0;
        0  0  1]
end

function analyticHcoeffs(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    Sym[1  0  0;
        0  1  0;
        0  0  1]
end


macro zpinch_equilibrium(R₀, B₀)
    generate_equilibrium_code(ThetaPinch(R₀, B₀); output=false)
end
