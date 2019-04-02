"""
θ-pinch equilibrium in (x,y,z) coordinates.

Parameters:
    B₀: B-field at magnetic axis
"""
struct ThetaPinch{T <: Number} <: CartesianEquilibrium
    name::String
    B₀::T

    function ThetaPinch{T}(B₀::T) where T <: Number
        new("ThetaPinchEquilibrium", B₀)
    end
end

ThetaPinch(B₀::T=1.0) where T <: Number = ThetaPinch{T}(B₀)


function Base.show(io::IO, equ::ThetaPinch)
    print(io, "θ-Pinch Equilibrium in (x,y,z) Coordinates with\n")
    print(io, "  B₀ = ", equ.B₀)
end


function R(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    r(x,equ)
end

function r(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    sqrt(r²(x,equ))
end

function r²(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    X(x,equ)^2 + Y(x,equ)^2
end

function θ(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    atan2(Y(x,equ), X(x,equ))
end

function ϕ(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    θ(x,equ)
end


function A₁(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    - equ.B₀ * Y(x,equ) / 2
end

function A₂(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    + equ.B₀ * X(x,equ) / 2
end

function A₃(x::AbstractArray{T,1}, equ::ThetaPinch) where {T <: Number}
    zero(eltype(x))
end


macro zpinch_equilibrium(R₀, B₀)
    generate_equilibrium_code(ThetaPinch(R₀, B₀); output=false)
end
