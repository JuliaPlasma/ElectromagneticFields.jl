"""
Axisymmetric tokamak equilibrium in (x,y,z) coordinates.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    q:  safety factor
"""
struct AxisymmetricTokamakCartesian{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q::T

    function AxisymmetricTokamakCartesian{T}(R₀::T, B₀::T, q::T) where T <: Number
        new("AxisymmetricTokamakCartesianEquilibrium", R₀, B₀, q)
    end
end

AxisymmetricTokamakCartesian(R₀::T=1.0, B₀::T=1.0, q::T=2.0) where T <: Number = AxisymmetricTokamakCartesian{T}(R₀, B₀, q)


function Base.show(io::IO, equ::AxisymmetricTokamakCartesian)
    print(io, "Axisymmetric Tokamak Equilibrium in (x,y,z) Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q  = ", equ.q)
end


function r²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    (R(x, equ) - equ.R₀)^2 + Z(x, equ)^2
end

function r(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    sqrt(r²(x, equ))
end

function θ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    atan2(Z(x, equ), R(x, equ) - equ.R₀)
end

function R(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    sqrt(x[1]^2 + x[2]^2)
end

function Z(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    x[3]
end

function ϕ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    atan2(x[2], x[1])
end


function analyticA₁(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    + 0.5 * equ.B₀ * equ.R₀ * Z(x, equ) / R(x, equ)
end

function analyticA₂(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    - 0.5 * equ.B₀ * equ.R₀ * log(R(x, equ) / equ.R₀)
end

function analyticA₃(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    - 0.5 * equ.B₀ * r²(x, equ) / (equ.q * R(x, equ))
end

function analyticMetric(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    [1  0  0;
     0  1  0;
     0  0  1]
end

function analyticHcoeffs(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCartesian) where {T <: Number}
    [1  0  0;
     0  1  0;
     0  0  1]
end


macro axisymmetric_tokamak_equilibrium_cartesian(R₀, B₀, q)
    generate_equilibrium_code(AxisymmetricTokamakCartesian(R₀, B₀, q); output=false)
end
