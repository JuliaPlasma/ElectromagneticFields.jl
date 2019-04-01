
abstract type AbstractSolovevEquilibrium <: AnalyticEquilibrium end



function X(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    R(x,equ) * cos(ϕ(x,equ))
end

function Y(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    R(x,equ) * sin(ϕ(x,equ))
end

function Z(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    x[2] * equ.R₀
end

function R(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    x[1] * equ.R₀
end

function r(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    sqrt(r²(x, equ))
end

function r²(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    (R(x,equ) - equ.R₀)^2 + Z(x,equ)^2
end

function θ(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    atan2(Z(x,equ), R(x,equ) - equ.R₀)
end

function ϕ(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    x[3]
end


function J(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    R(x,equ)
end


function periodicity(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    p = zero(x)
    p[3] = 2π
    return p
end


function A₁(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    + equ.B₀ * equ.R₀ * x[2] / x[1] / 2
end

function A₂(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    - equ.B₀ * equ.R₀ * log(x[1]) / 2
end


function g₁₁(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    one(T)
end

function g₂₂(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    one(T)
end

function g₃₃(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    R(x, equ)^2
end
