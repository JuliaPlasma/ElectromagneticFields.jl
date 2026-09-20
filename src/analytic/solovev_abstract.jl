
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

J(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x, equ) * equ.R₀^2
# (R/R₀, Z/R₀, ϕ) is the same left-handed ordering as the cylindrical chart.
orientation(::AbstractSolovevEquilibrium) = -1

function A₁(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    +equ.B₀ * equ.R₀ * x[2] / x[1] / 2
end
function A₂(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    -equ.B₀ * equ.R₀ * NaNMath.log(x[1]) / 2
end

x¹(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = X(ξ, equ)
x²(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = Y(ξ, equ)
x³(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = Z(ξ, equ)

function ξ¹(x::AbstractVector, equ::AbstractSolovevEquilibrium)
    sqrt(x[1]^2 + x[2]^2) / equ.R₀
end
ξ²(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[3] / equ.R₀
ξ³(x::AbstractVector, equ::AbstractSolovevEquilibrium) = atan(x[2], x[1])

g₁₁(x::AbstractVector, equ::AbstractSolovevEquilibrium) = equ.R₀^2
g₂₂(x::AbstractVector, equ::AbstractSolovevEquilibrium) = equ.R₀^2
g₃₃(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x, equ)^2

function get_functions(::AbstractSolovevEquilibrium)
    (X = X, Y = Y, Z = Z, R = R, r = r, θ = θ, ϕ = ϕ, r² = r²)
end
# `c` is deliberately not excluded here, although it is derived from the other six rather than
# chosen: `A₃` reads it, so it has to reach the generated code as a parameter like the rest. The
# default — every field but `name` — is therefore what this type wants, and no `get_parameters`
# method is defined for it.

minx³(ξ::AbstractVector{T}, equ::AbstractSolovevEquilibrium) where {T} = T(0)
maxx³(ξ::AbstractVector{T}, equ::AbstractSolovevEquilibrium) where {T} = T(2π)

# (R/R₀, Z/R₀, ϕ) is the cylindrical chart, so the toroidal angle alone is periodic.
GeometricBase.periodic(::AbstractSolovevEquilibrium) = SVector(false, false, true)
