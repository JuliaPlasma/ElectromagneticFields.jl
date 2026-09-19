
const DEFAULT_DIPOLE_B₀ = 1000.0

@doc raw"""
Dipole magnetic field in (x,y,z) coordinates
Based on Xinjie Li, Ruili Zhang, and Jian Liu, Symplectic Runge-Kutta methods for the guiding
    center dynamics.

The covariant components of the vector potential are given by
```math
A (x,y,z) = \frac{B₀}{r^3} \big( y , \, -x , \, 0 \big)^T ,
```
resulting in the magnetic field with covariant components
```math
B (x,y,z) = - \frac{B₀}{r^5} \big( 3xz, \, 3yz, \, 2z^2 - x^2 - y^2 \big)^T .
```
"""
struct DipoleField{T <: Number} <: CartesianEquilibrium
    name::String
    B₀::T

    function DipoleField{T}(B₀::T) where {T <: Number}
        new("DipoleField", B₀)
    end
end

DipoleField(B₀::T = DEFAULT_DIPOLE_B₀) where {T} = DipoleField{T}(B₀)

function Base.show(io::IO, equ::DipoleField)
    print(io, "Dipole Field in (x,y,z) Coordinates")
end

function A₁(x::AbstractVector, equ::DipoleField)
    +equ.B₀ * Y(x, equ) / sqrt(X(x, equ)^2 + Y(x, equ)^2 + Z(x, equ)^2)^3
end
function A₂(x::AbstractVector, equ::DipoleField)
    -equ.B₀ * X(x, equ) / sqrt(X(x, equ)^2 + Y(x, equ)^2 + Z(x, equ)^2)^3
end
A₃(x::AbstractVector, equ::DipoleField) = zero(eltype(x))

get_functions(::DipoleField) = (X = X, Y = Y, Z = Z)
