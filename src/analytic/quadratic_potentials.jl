@doc raw"""
Electromagnetic field with quadratic potentials in (x,y,z) coordinates
Based on Xinjie Li, Ruili Zhang, and Jian Liu, Symplectic Runge-Kutta methods for the guiding
    center dynamics.

The covariant components of the vector potential are given by
```math
A (x,y,z) = \bigg( -50y , \, 50x , \, \frac{x^2 + y^2}{2} \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (x,y,z) = \big( y, \, -x, \, 100 \big)^T ,
```
and electrostatic potential
```math
\phi (x,y,z) = \frac{1}{2} \big( x^2 + y^2 + z^2 \big) .
```
"""
module QuadraticPotentials


import ..ElectromagneticFields
import ..ElectromagneticFields: CartesianEquilibrium, code, code_arguments
import ..AnalyticCartesianField: X, Y, Z

export QuadraticPotentialsField

const DEFAULT_Bz = 100.0

struct QuadraticPotentialsField{T<:Number} <: CartesianEquilibrium
    name::String
    Bz::T

    function QuadraticPotentialsField{T}(Bz::T) where {T<:Number}
        new("QuadraticPotentialsField", Bz)
    end
end

QuadraticPotentialsField(Bz::T=DEFAULT_Bz) where {T} = QuadraticPotentialsField{T}(Bz)

function init(args...)
    QuadraticPotentialsField(args...)
end

macro code(args...)
    parameters, options = code_arguments(args)
    code(init(parameters...); escape=true, options...)
end

function Base.show(io::IO, equ::QuadraticPotentialsField)
    print(io, "Electromagnetic field with quadratic potentials in (x,y,z) coordinates")
end


ElectromagneticFields.A₁(x::AbstractVector, equ::QuadraticPotentialsField) = -equ.Bz * Y(x, equ) / 2
ElectromagneticFields.A₂(x::AbstractVector, equ::QuadraticPotentialsField) = +equ.Bz * X(x, equ) / 2
ElectromagneticFields.A₃(x::AbstractVector, equ::QuadraticPotentialsField) = (X(x, equ)^2 + Y(x, equ)^2) / 2
ElectromagneticFields.φ(x::AbstractVector, equ::QuadraticPotentialsField) = (X(x, equ)^2 + Y(x, equ)^2 + Z(x, equ)^2) / 2

ElectromagneticFields.get_functions(::QuadraticPotentialsField) = (X=X, Y=Y, Z=Z)


end
