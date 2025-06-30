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

using RecipesBase
using LaTeXStrings

import ..ElectromagneticFields
import ..ElectromagneticFields: CartesianEquilibrium, code
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
    code(init(args...); escape=true)
end

function Base.show(io::IO, equ::QuadraticPotentialsField)
    print(io, "Electromagnetic field with quadratic potentials in (x,y,z) coordinates")
end


ElectromagneticFields.A₁(x::AbstractVector, equ::QuadraticPotentialsField) = -equ.Bz * Y(x, equ) / 2
ElectromagneticFields.A₂(x::AbstractVector, equ::QuadraticPotentialsField) = +equ.Bz * X(x, equ) / 2
ElectromagneticFields.A₃(x::AbstractVector, equ::QuadraticPotentialsField) = (X(x, equ)^2 + Y(x, equ)^2) / 2
ElectromagneticFields.φ(x::AbstractVector, equ::QuadraticPotentialsField) = (X(x, equ)^2 + Y(x, equ)^2 + Z(x, equ)^2) / 2

ElectromagneticFields.get_functions(::QuadraticPotentialsField) = (X=X, Y=Y, Z=Z)


@recipe function f(equ::QuadraticPotentialsField;
    nx=100, ny=100, levels=20, size=(1200, 400),
    xlims=(-1.0, +1.0),
    ylims=(-1.0, +1.0))

    xgrid = LinRange(xlims[1], xlims[2], nx)
    ygrid = LinRange(ylims[1], ylims[2], ny)
    pot1 = [ElectromagneticFields.A₁([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
    pot2 = [ElectromagneticFields.A₂([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
    pot3 = [ElectromagneticFields.A₃([xgrid[i], ygrid[j], 0.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]

    seriestype := :contour
    aspect_ratio := :equal
    layout := (1, 3)
    size := size
    xlims := xlims
    ylims := ylims
    levels := levels
    legend := :none

    @series begin
        subplot := 1
        title := L"A_x (x,y,0)"
        xguide := L"x"
        yguide := L"y"
        (xgrid, ygrid, pot1)
    end

    @series begin
        subplot := 2
        title := L"A_y (x,y,0)"
        xguide := L"x"
        yguide := L"y"
        (xgrid, ygrid, pot2)
    end

    @series begin
        subplot := 3
        title := L"A_z (x,y,0)"
        xguide := L"x"
        yguide := L"y"
        (xgrid, ygrid, pot3)
    end
end

end
