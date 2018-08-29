
abstract type AbstractSolovevEquilibrium <: AnalyticEquilibrium end


function analyticA₁(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    + equ.B₀ * x[2] / x[1] / 2
end

function analyticA₂(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    - equ.B₀ * equ.R₀ * log(x[1]) / 2
end


function analyticMetric(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    R = x[1] * equ.R₀
    [1  0  0;
     0  1  0;
     0  0  R^2]
end

function analyticHcoeffs(x::AbstractArray{T,1}, equ::AbstractSolovevEquilibrium) where {T <: Number}
    [1  0  0;
     0  1  0;
     0  0  x[1] * equ.R₀]
end
