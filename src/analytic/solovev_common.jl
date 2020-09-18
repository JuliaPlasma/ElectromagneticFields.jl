
abstract type AbstractSolovevEquilibrium <: AnalyticEquilibrium end


Z(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[2] * equ.R₀
R(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[1] * equ.R₀
ϕ(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[3]

X(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x,equ) * cos(ϕ(x,equ))
Y(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x,equ) * sin(ϕ(x,equ))
θ(x::AbstractVector, equ::AbstractSolovevEquilibrium) = atan(Z(x,equ), R(x,equ) - equ.R₀)

r²(x::AbstractVector, equ::AbstractSolovevEquilibrium) = (R(x,equ) - equ.R₀)^2 + Z(x,equ)^2
r(x::AbstractVector, equ::AbstractSolovevEquilibrium) = sqrt(r²(x, equ))

J(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x,equ)

A₁(x::AbstractVector, equ::AbstractSolovevEquilibrium) = + equ.B₀ * equ.R₀ * x[2] / x[1] / 2
A₂(x::AbstractVector, equ::AbstractSolovevEquilibrium) = - equ.B₀ * equ.R₀ * log(x[1]) / 2

x¹(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = X(ξ,equ)
x²(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = Y(ξ,equ)
x³(ξ::AbstractVector, equ::AbstractSolovevEquilibrium) = Z(ξ,equ)

ξ¹(x::AbstractVector, equ::AbstractSolovevEquilibrium) = sqrt(x[1]^2 + x[2]^2) / equ.R₀
ξ²(x::AbstractVector, equ::AbstractSolovevEquilibrium) = x[3] / equ.R₀
ξ³(x::AbstractVector, equ::AbstractSolovevEquilibrium) = atan(x[2], x[1])

g₁₁(x::AbstractVector, equ::AbstractSolovevEquilibrium) = one(eltype(x))
g₂₂(x::AbstractVector, equ::AbstractSolovevEquilibrium) = one(eltype(x))
g₃₃(x::AbstractVector, equ::AbstractSolovevEquilibrium) = R(x, equ)^2

get_functions(::AbstractSolovevEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)

function periodicity(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    p = zero(x)
    p[3] = 2π
    return p
end
