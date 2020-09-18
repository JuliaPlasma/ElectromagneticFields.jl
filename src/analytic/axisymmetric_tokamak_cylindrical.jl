@doc raw"""
Axisymmetric tokamak equilibrium in (R,Z,ϕ) coordinates with covariant
components of the vector potential given by
```math
A (R, Z, \phi) = \frac{B_0}{2} \, \bigg( R_0 \, \frac{Z}{R} , \, - R_0 \, \ln \bigg( \frac{R}{R_0} \bigg) , \, - \frac{r^2}{q_0} \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (R, Z, \phi) = \frac{B_0}{q_0} \, \bigg( - \frac{Z}{R} , \, \frac{R - R_0}{R} , \, - q_0 R_0 \bigg)^T ,
```
where $r = \sqrt{ (R - R_0)^2 + Z^2 }$.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `q₀`: safety factor at magnetic axis
"""
struct AxisymmetricTokamakCylindrical{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakCylindrical{T}(R₀::T, B₀::T, q₀::T) where T <: Number
        new("AxisymmetricTokamakCylindricalEquilibrium", R₀, B₀, q₀)
    end
end

AxisymmetricTokamakCylindrical(R₀::T=1.0, B₀::T=1.0, q₀::T=2.0) where T <: Number = AxisymmetricTokamakCylindrical{T}(R₀, B₀, q₀)


function Base.show(io::IO, equ::AxisymmetricTokamakCylindrical)
    print(io, "Axisymmetric Tokamak Equilibrium in (R,Z,ϕ) Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q₀ = ", equ.q₀)
end


R(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = x[1]
Z(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = x[2]
ϕ(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = x[3]
r²(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = (R(x,equ) - equ.R₀)^2 + Z(x,equ)^2
r(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = sqrt(r²(x, equ))
X(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = R(x,equ) * cos(ϕ(x,equ))
Y(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = R(x,equ) * sin(ϕ(x,equ))
θ(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = atan(Z(x,equ), R(x,equ) - equ.R₀)

J(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = R(x,equ)

A₁(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = + equ.B₀ * equ.R₀ * Z(x,equ) / R(x,equ) / 2
A₂(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = - equ.B₀ * equ.R₀ * log(R(x,equ) / equ.R₀) / 2
A₃(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = - equ.B₀ * r²(x,equ) / equ.q₀ / 2

g₁₁(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = one(eltype(x))
g₂₂(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = one(eltype(x))
g₃₃(x::AbstractVector, equ::AxisymmetricTokamakCylindrical) = R(x, equ)^2


get_functions(::AxisymmetricTokamakCylindrical) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)

function periodicity(x::AbstractVector, ::AxisymmetricTokamakCylindrical)
    p = zero(x)
    p[3] = 2π
    return p
end


macro axisymmetric_tokamak_equilibrium_cylindrical(R₀, B₀, q₀)
    generate_equilibrium_code(AxisymmetricTokamakCylindrical(R₀, B₀, q₀); output=false)
end
