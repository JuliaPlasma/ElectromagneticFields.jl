
abstract type CartesianEquilibrium <: AnalyticEquilibrium end
abstract type CartesianPerturbation <: AnalyticPerturbation end

CartesianEquPert = Union{CartesianEquilibrium, CartesianPerturbation}


function X(x::AbstractArray{T,1}, equ::CartesianEquPert) where {T <: Number}
    x[1]
end

function Y(x::AbstractArray{T,1}, equ::CartesianEquPert) where {T <: Number}
    x[2]
end

function Z(x::AbstractArray{T,1}, equ::CartesianEquPert) where {T <: Number}
    x[3]
end

function J(x::AbstractArray{T,1}, equ::CartesianEquPert) where {T <: Number}
    one(T)
end

function g₁₁(x::AbstractArray{T,1}, equ::CartesianEquPert) where {T <: Number}
    one(T)
end

function g₂₂(x::AbstractArray{T,1}, equ::CartesianEquPert) where {T <: Number}
    one(T)
end

function g₃₃(x::AbstractArray{T,1}, equ::CartesianEquPert) where {T <: Number}
    one(T)
end
