
"""
Axisymmetric tokamak equilibrium in (r,θ,ϕ) coordinates.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    q:  safety factor
"""
struct AxisymmetricTokamakToroidal{T <: Number} <: AnalyticEquilibrium
    const name::String = "AxisymmetricTokamakEquilibriumToroidal"
    R₀::T
    B₀::T
    q::T

    function AxisymmetricTokamakToroidal{T}(R₀::T, B₀::T, q::T) where T <: Number
        new(R₀, B₀, q)
    end
end

AxisymmetricTokamakToroidal(R₀::T, B₀::T, q::T) where T <: Number = AxisymmetricTokamakToroidal{T}(R₀, B₀, q)


function Base.show(io::IO, equ::AxisymmetricTokamakToroidal)
    print(io, "Axisymmetric Tokamak Equilibrium in Toroidal Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q  = ", equ.q)
end


function r²(x::Vector, equ::AxisymmetricTokamakToroidal)
    x[1]^2
end

function r(x::Vector, equ::AxisymmetricTokamakToroidal)
    x[1]
end

function θ(x::Vector, equ::AxisymmetricTokamakToroidal)
    x[2]
end

function R(x::Vector, equ::AxisymmetricTokamakToroidal)
    equ.R₀ + x[1] * cos(x[2])
end

function Z(x::Vector, equ::AxisymmetricTokamakToroidal)
    x[1] * sin(x[2])
end

function ϕ(x::Vector, equ::AxisymmetricTokamakToroidal)
    x[3]
end


function analyticA₁(x::Vector, equ::AxisymmetricTokamakToroidal)
    zero(x)
end

function analyticA₂(x::Vector, equ::AxisymmetricTokamakToroidal)
    equ.B₀ * equ.R₀ / cos(x[2])^2 * ( x[1] * cos(x[2]) - equ.R₀ * log(one(x) + x[1] * cos(x[2]) / equ.R₀) )
end

function analyticA₃(x::Vector, equ::AxisymmetricTokamakToroidal)
    - 0.5 * equ.B₀ * r²(x, equ) / equ.q
end

function analyticMetric(x::Vector, equ::AxisymmetricTokamakToroidal)
    [1   0             0;
     0   r(x, equ)^2   0;
     0   0             R(x, equ)^2]
end


macro axisymmetric_tokamak_equilibrium_toroidal(R₀, B₀, q)
    generate_equilibrium_code(AxisymmetricTokamakToroidal(R₀, B₀, q); output=false)
end
