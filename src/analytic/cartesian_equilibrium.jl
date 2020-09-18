
abstract type CartesianEquilibrium <: AnalyticEquilibrium end
abstract type CartesianPerturbation <: AnalyticPerturbation end

CartesianField = Union{CartesianEquilibrium, CartesianPerturbation}


X(x::AbstractVector, ::CartesianField) = x[1]
Y(x::AbstractVector, ::CartesianField) = x[2]
Z(x::AbstractVector, ::CartesianField) = x[3]

J(x::AbstractVector, ::CartesianField) = one(eltype(x))

x¹(ξ::AbstractVector, ::CartesianField) = ξ[1]
x²(ξ::AbstractVector, ::CartesianField) = ξ[2]
x³(ξ::AbstractVector, ::CartesianField) = ξ[3]

ξ¹(x::AbstractVector, ::CartesianField) = x[1]
ξ²(x::AbstractVector, ::CartesianField) = x[2]
ξ³(x::AbstractVector, ::CartesianField) = x[3]

g₁₁(x::AbstractVector, ::CartesianField) = one(eltype(x))
g₂₂(x::AbstractVector, ::CartesianField) = one(eltype(x))
g₃₃(x::AbstractVector, ::CartesianField) = one(eltype(x))
