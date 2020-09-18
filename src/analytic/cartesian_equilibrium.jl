
abstract type CartesianEquilibrium <: AnalyticEquilibrium end
abstract type CartesianPerturbation <: AnalyticPerturbation end

CartesianField = Union{CartesianEquilibrium, CartesianPerturbation}


X(x::AbstractVector, ::CartesianField) = x[1]
Y(x::AbstractVector, ::CartesianField) = x[2]
Z(x::AbstractVector, ::CartesianField) = x[3]

J(x::AbstractVector, ::CartesianField) = one(eltype(x))

g₁₁(x::AbstractVector, ::CartesianField) = one(eltype(x))
g₂₂(x::AbstractVector, ::CartesianField) = one(eltype(x))
g₃₃(x::AbstractVector, ::CartesianField) = one(eltype(x))
