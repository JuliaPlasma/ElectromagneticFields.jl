"""
Axisymmetric tokamak equilibrium in (r,θ,ϕ) coordinates.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    q:  safety factor
"""
struct AxisymmetricTokamakToroidal{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q::T

    function AxisymmetricTokamakToroidal{T}(R₀::T, B₀::T, q::T) where T <: Number
        new("AxisymmetricTokamakEquilibriumToroidal", R₀, B₀, q)
    end
end

AxisymmetricTokamakToroidal(R₀::T=1.0, B₀::T=1.0, q::T=2.0) where T <: Number = AxisymmetricTokamakToroidal{T}(R₀, B₀, q)


function Base.show(io::IO, equ::AxisymmetricTokamakToroidal)
    print(io, "Axisymmetric Tokamak Equilibrium in Toroidal Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q  = ", equ.q)
end


function r(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    x[1]
end

function θ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    x[2]
end

function R(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    equ.R₀ + x[1] * cos(x[2])
end

function Z(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    x[1] * sin(x[2])
end

function ϕ(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    x[3]
end

function r²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    x[1]^2
end

function R²(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    expand( (equ.R₀ + x[1] * cos(x[2]))^2 )
end


function analyticA₁(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    zero(T)
end

function analyticA₂(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    equ.B₀ * equ.R₀ / cos(x[2])^2 * ( x[1] * cos(x[2]) - equ.R₀ * log(x[1] * cos(x[2]) / equ.R₀ + 1) )
end

function analyticA₃(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    - 0.5 * equ.B₀ * r²(x, equ) / equ.q
end

function analyticMetric(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    [1   0            0;
     0   r²(x, equ)   0;
     0   0            R²(x, equ)]
end

function analyticHcoeffs(x::AbstractArray{T,1}, equ::AxisymmetricTokamakToroidal) where {T <: Number}
    [1   0           0;
     0   r(x, equ)   0;
     0   0           R(x, equ)]
end


macro axisymmetric_tokamak_equilibrium_toroidal(R₀, B₀, q)
    generate_equilibrium_code(AxisymmetricTokamakToroidal(R₀, B₀, q); output=false)
end
