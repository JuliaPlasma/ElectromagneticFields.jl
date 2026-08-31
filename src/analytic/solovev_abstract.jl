module SolovevAbstract

import NaNMath: log

import ..ElectromagneticFields
import ..ElectromagneticFields: AnalyticEquilibrium, AnalyticPerturbation

export AbstractSolovevEquilibrium

abstract type AbstractSolovevEquilibrium <: AnalyticEquilibrium end

R(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[1] * equ.R₀
Z(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[2] * equ.R₀
ϕ(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[3]

X(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x, equ) * cos(ϕ(x, equ))
Y(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x, equ) * sin(ϕ(x, equ))
θ(x::AbstractVector, equ::AbstractSolovevEquilibrium) = atan(Z(x, equ), R(x, equ) - equ.R₀)

function r²(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    (R(x, equ) - equ.R₀)^2 + Z(x, equ)^2
end
r(x::AbstractVector, equ::AbstractSolovevEquilibrium) = sqrt(r²(x, equ))

function ElectromagneticFields.J(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    R(x, equ) * equ.R₀^2
end
# (R/R₀, Z/R₀, ϕ) is the same left-handed ordering as the cylindrical chart.
ElectromagneticFields.orientation(::AbstractSolovevEquilibrium) = -1

function ElectromagneticFields.A₁(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    +equ.B₀ * equ.R₀ * x[2] / x[1] / 2
end
function ElectromagneticFields.A₂(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    -equ.B₀ * equ.R₀ * log(x[1]) / 2
end

ElectromagneticFields.x¹(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = X(ξ, equ)
ElectromagneticFields.x²(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = Y(ξ, equ)
ElectromagneticFields.x³(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = Z(ξ, equ)

function ElectromagneticFields.ξ¹(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    sqrt(x[1]^2 + x[2]^2) / equ.R₀
end
ElectromagneticFields.ξ²(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[3] / equ.R₀
function ElectromagneticFields.ξ³(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    atan(x[2], x[1])
end

ElectromagneticFields.g₁₁(x::AbstractVector, equ::AbstractSolovevEquilibrium) = equ.R₀^2
ElectromagneticFields.g₂₂(x::AbstractVector, equ::AbstractSolovevEquilibrium) = equ.R₀^2
ElectromagneticFields.g₃₃(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x, equ)^2

function ElectromagneticFields.get_functions(::AbstractSolovevEquilibrium)
    (X = X, Y = Y, Z = Z, R = R, r = r, θ = θ, ϕ = ϕ, r² = r²)
end
function ElectromagneticFields.get_parameters(::AbstractSolovevEquilibrium)
    (:R₀, :B₀, :ϵ, :κ, :δ, :α)
end

function ElectromagneticFields.minx³(ξ::AbstractVector{T}, equ::AbstractSolovevEquilibrium) where {T}
    T(0)
end
function ElectromagneticFields.maxx³(ξ::AbstractVector{T}, equ::AbstractSolovevEquilibrium) where {T}
    T(2π)
end

end
