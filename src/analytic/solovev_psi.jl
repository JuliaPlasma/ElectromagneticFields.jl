
function ψ₀(x::AbstractVector{T}, a) where {T <: Number}
    x[1]^4 / 8 + a * (x[1]^2 * log(x[1]) / 2 - x[1]^4 / 8 )
end

function ψ₁(x::AbstractVector{T}) where {T <: Number}
    one(T)
end

function ψ₂(x::AbstractVector{T}) where {T <: Number}
    x[1]^2
end

function ψ₃(x::AbstractVector{T}) where {T <: Number}
    x[2]^2 - x[1]^2 * log(x[1])
end

function ψ₄(x::AbstractVector{T}) where {T <: Number}
    x[1]^4 - 4 * x[1]^2 * x[2]^2
end

function ψ₅(x::AbstractVector{T}) where {T <: Number}
    2 * x[2]^4 - 9 * x[2]^2 * x[1]^2 + 3 * x[1]^4 * log(x[1]) - 12 * x[1]^2 * x[2]^2 * log(x[1])
end

function ψ₆(x::AbstractVector{T}) where {T <: Number}
    x[1]^6 - 12 * x[1]^4 * x[2]^2 + 8 * x[1]^2 * x[2]^4
end

function ψ₇(x::AbstractVector{T}) where {T <: Number}
    8 * x[2]^6 - 140 * x[2]^4 * x[1]^2 + 75 * x[2]^2 * x[1]^4 - 15 * x[1]^6 * log(x[1]) +
                 180 * x[1]^4 * x[2]^2 * log(x[1]) - 120 * x[1]^2 * x[2]^4 * log(x[1])
end

function ψ₈(x::AbstractVector{T}) where {T <: Number}
    x[2]
end

function ψ₉(x::AbstractVector{T}) where {T <: Number}
    x[2] * x[1]^2
end

function ψ₁₀(x::AbstractVector{T}) where {T <: Number}
    x[2]^3 - 3 * x[2] * x[1]^2 * log(x[1])
end

function ψ₁₁(x::AbstractVector{T}) where {T <: Number}
    3 * x[2] * x[1]^4 - 4 * x[2]^3 * x[1]^2
end

function ψ₁₂(x::AbstractVector{T}) where {T <: Number}
    8 * x[2]^5 - 45 * x[2] * x[1]^4 - 80 * x[2]^3 * x[1]^2 * log(x[1]) + 60 * x[2] * x[1]^4 * log(x[1])
end
 