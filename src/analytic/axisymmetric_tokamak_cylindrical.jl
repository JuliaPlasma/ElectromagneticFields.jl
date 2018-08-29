"""
Axisymmetric tokamak equilibrium in (R,Z,ϕ) coordinates.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    q:  safety factor
"""
struct AxisymmetricTokamakCylindrical{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q::T

    function AxisymmetricTokamakCylindrical{T}(R₀::T, B₀::T, q::T) where T <: Number
        new("AxisymmetricTokamakCylindricalEquilibrium", R₀, B₀, q)
    end
end

AxisymmetricTokamakCylindrical(R₀::T=1.0, B₀::T=1.0, q::T=2.0) where T <: Number = AxisymmetricTokamakCylindrical{T}(R₀, B₀, q)


function Base.show(io::IO, equ::AxisymmetricTokamakCylindrical)
    print(io, "Axisymmetric Tokamak Equilibrium in (R,Z,ϕ) Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q  = ", equ.q)
end


function r²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    (x[1] - equ.R₀)^2 + x[2]^2
end

function r(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    sqrt(r²(x, equ))
end

function θ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    atan2(x[2], x[1] - equ.R₀)
end

function R(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    x[1]
end

function Z(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    x[2]
end

function ϕ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    x[3]
end


function analyticA₁(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    + 0.5 * equ.B₀ * equ.R₀ * x[2] / x[1]
end

function analyticA₂(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    - 0.5 * equ.B₀ * equ.R₀ * log(x[1] / equ.R₀)
end

function analyticA₃(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    - 0.5 * equ.B₀ * r²(x, equ) / (equ.q * x[1])
end

function analyticMetric(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    [1  0  0;
     0  1  0;
     0  0  x[1]^2]
end

function analyticHcoeffs(x::AbstractArray{T,1}, equ::AxisymmetricTokamakCylindrical) where {T <: Number}
    [1  0  0;
     0  1  0;
     0  0  x[1]]
end


macro axisymmetric_tokamak_equilibrium_cylindrical(R₀, B₀, q)
    generate_equilibrium_code(AxisymmetricTokamakCylindrical(R₀, B₀, q); output=false)
end
