
"""
Axisymmetric tokamak equilibrium in (R,Z,ϕ) coordinates.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    q:  safety factor
"""
struct AxisymmetricTokamak{T <: Number} <: AnalyticEquilibrium
    const name::String = "AxisymmetricTokamakEquilibrium"
    R₀::T
    B₀::T
    q::T

    function AxisymmetricTokamak{T}(R₀::T, B₀::T, q::T) where T <: Number
        new(R₀, B₀, q)
    end
end

AxisymmetricTokamak(R₀::T, B₀::T, q::T) where T <: Number = AxisymmetricTokamak{T}(R₀, B₀, q)


function Base.show(io::IO, equ::AxisymmetricTokamak)
    print(io, "Axisymmetric Tokamak Equilibrium with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q  = ", equ.q)
end


function r²(x::Vector, equ::AxisymmetricTokamak)
    (x[1] - equ.R₀)^2 + x[2]^2
end

function r(x::Vector, equ::AxisymmetricTokamak)
    sqrt(r²(x, equ))
end

function θ(x::Vector, equ::AxisymmetricTokamak)
    atan2(x[2], x[1] - equ.R₀)
end

function R(x::Vector, equ::AxisymmetricTokamak)
    x[1]
end

function Z(x::Vector, equ::AxisymmetricTokamak)
    x[2]
end

function ϕ(x::Vector, equ::AxisymmetricTokamak)
    x[3]
end


function analyticA₁(x::Vector, equ::AxisymmetricTokamak)
    + 0.5 * equ.B₀ * equ.R₀ * x[2] / x[1]
end

function analyticA₂(x::Vector, equ::AxisymmetricTokamak)
    - 0.5 * equ.B₀ * equ.R₀ * log(x[1] / equ.R₀)
end

function analyticA₃(x::Vector, equ::AxisymmetricTokamak)
    - 0.5 * equ.B₀ * r²(x, equ) / (equ.q * x[1])
end

function analyticMetric(x::Vector, equ::AxisymmetricTokamak)
    [1  0  0;
     0  1  0;
     0  0  x[1]^2]
end


macro axisymmetric_tokamak_equilibrium(R₀, B₀, q)
    generate_equilibrium_code(AxisymmetricTokamak(R₀, B₀, q); output=false)
end
