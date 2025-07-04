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
module Dipole

using RecipesBase
using LaTeXStrings

import ..ElectromagneticFields
import ..ElectromagneticFields: CartesianEquilibrium, code
import ..AnalyticCartesianField: X, Y, Z

export DipoleField

const DEFAULT_B₀ = 1000.0

struct DipoleField{T<:Number} <: CartesianEquilibrium
    name::String
    B₀::T

    function DipoleField{T}(B₀::T) where {T<:Number}
        new("DipoleField", B₀)
    end
end

DipoleField(B₀::T=DEFAULT_B₀) where {T} = DipoleField{T}(B₀)

function init(args...)
    DipoleField(args...)
end

macro code(args...)
    code(init(args...); escape=true)
end

function Base.show(io::IO, equ::DipoleField)
    print(io, "Dipole Field in (x,y,z) Coordinates")
end


ElectromagneticFields.A₁(x::AbstractVector, equ::DipoleField) = +equ.B₀ * Y(x, equ) / sqrt(X(x, equ)^2 + Y(x, equ)^2 + Z(x, equ)^2)^3
ElectromagneticFields.A₂(x::AbstractVector, equ::DipoleField) = -equ.B₀ * X(x, equ) / sqrt(X(x, equ)^2 + Y(x, equ)^2 + Z(x, equ)^2)^3
ElectromagneticFields.A₃(x::AbstractVector, equ::DipoleField) = zero(eltype(x))

ElectromagneticFields.get_functions(::DipoleField) = (X=X, Y=Y, Z=Z)


@recipe function f(equ::DipoleField;
    nx=100, ny=100, levels=20, size=(800, 400),
    xlims=(-1.0, +1.0),
    ylims=(-1.0, +1.0))

    xgrid = LinRange(xlims[1], xlims[2], nx)
    ygrid = LinRange(ylims[1], ylims[2], ny)
    pot1 = [ElectromagneticFields.A₁([xgrid[i], ygrid[j], 1.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]
    pot2 = [ElectromagneticFields.A₂([xgrid[i], ygrid[j], 1.0], equ) for i in eachindex(xgrid), j in eachindex(ygrid)]

    seriestype := :contour
    aspect_ratio := :equal
    layout := (1, 2)
    size := size
    xlims := xlims
    ylims := ylims
    levels := levels
    legend := :none

    @series begin
        subplot := 1
        title := L"A_x (x,y,1)"
        xguide := L"x"
        yguide := L"y"
        (xgrid, ygrid, pot1)
    end

    @series begin
        subplot := 2
        title := L"A_y (x,y,1)"
        xguide := L"x"
        yguide := L"y"
        (xgrid, ygrid, pot2)
    end
end

end
