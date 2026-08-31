@doc raw"""
Symmetric Solov'ev equilibrium in cartesian (x,y,z) coordinates.
Based on McCarthy, Physics of Plasmas 6, 3554, 1999.

The covariant components of the vector potential are given by
```math
A (x, y) = \frac{B_0}{2} \, \bigg( 0 , \, 0 , \, - \frac{\alpha}{4} (R_0 + x)^4 - \beta y^2 \bigg)^T ,
```

Parameters:
 * `R₀`: major radius, which places the magnetic axis at `x = -R₀`
 * `B₀`: B-field at magnetic axis
 * `α`, `β`: free constants
"""
module SolovevSymmetric

import ..ElectromagneticFields
import ..ElectromagneticFields: CartesianEquilibrium, code, code_arguments
import ..AnalyticCartesianField: X, Y, Z

export SolovevSymmetricEquilibrium

const DEFAULT_R₀ = 1.0
const DEFAULT_B₀ = 1.0
const DEFAULT_α = 2.0
const DEFAULT_β = 0.5

struct SolovevSymmetricEquilibrium{T <: Number} <: CartesianEquilibrium
    name::String
    R₀::T
    B₀::T
    α::T
    β::T

    function SolovevSymmetricEquilibrium{T}(R₀::T, B₀::T, α::T, β::T) where {T <: Number}
        new("QuadraticSolovevEquilibrium", R₀, B₀, α, β)
    end
end

function SolovevSymmetricEquilibrium(R₀::T = DEFAULT_R₀, B₀::T = DEFAULT_B₀,
        α::T = DEFAULT_α, β::T = DEFAULT_β) where {T <: Number}
    SolovevSymmetricEquilibrium{T}(R₀, B₀, α, β)
end

function init(R₀ = DEFAULT_R₀, B₀ = DEFAULT_B₀, α = DEFAULT_α, β = DEFAULT_β)
    SolovevSymmetricEquilibrium(R₀, B₀, α, β)
end

macro code(args...)
    parameters, options = code_arguments(args)
    code(init(parameters...); escape = true, options...)
end

function Base.show(io::IO, equ::SolovevSymmetricEquilibrium)
    print(io, "Quadratic Solovev Equilibrium with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  α  = ", equ.α, "\n")
    print(io, "  β  = ", equ.β)
end

ElectromagneticFields.A₁(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = zero(eltype(x))
ElectromagneticFields.A₂(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = zero(eltype(x))
ElectromagneticFields.A₃(x::AbstractVector, equ::SolovevSymmetricEquilibrium) = - equ.B₀ *
                                                                                (equ.α *
                                                                                 (equ.R₀ +
                                                                                  X(x, equ))^4 /
                                                                                 4 +
                                                                                 equ.β *
                                                                                 Y(x, equ)^2) /
                                                                                2

ElectromagneticFields.get_functions(::SolovevSymmetricEquilibrium) = (X = X, Y = Y, Z = Z)

end
