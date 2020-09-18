using SymEngine

@doc raw"""
Axisymmetric tokamak equilibrium in (r,θ,ϕ) coordinates with covariant
components of the vector potential given by
```math
A (r, \theta, \phi) = B_0 \, \bigg( 0 , \, \frac{r R_0}{\cos (\theta)} - \bigg( \frac{R_0}{\cos (\theta)} \bigg)^2 \, \ln \bigg( \frac{R}{R_0} \bigg) , \, - \frac{r^2}{2 q_0} \bigg)^T ,
```
resulting in the magnetic field with covariant components
```math
B (r, \theta, \phi) = \frac{B_0}{q_0} \, \bigg( 0 , \, \frac{r^2}{R}, \, q_0 R_0 \bigg)^T ,
```
where $R = R_0 + r \cos \theta$.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `q₀`: safety factor at magnetic axis
"""
struct AxisymmetricTokamakToroidal{T <: Number} <: AnalyticEquilibrium
    name::String
    R₀::T
    B₀::T
    q₀::T

    function AxisymmetricTokamakToroidal{T}(R₀::T, B₀::T, q₀::T) where T <: Number
        new("AxisymmetricTokamakEquilibriumToroidal", R₀, B₀, q₀)
    end
end

AxisymmetricTokamakToroidal(R₀::T=1.0, B₀::T=1.0, q₀::T=2.0) where T <: Number = AxisymmetricTokamakToroidal{T}(R₀, B₀, q₀)


function Base.show(io::IO, equ::AxisymmetricTokamakToroidal)
    print(io, "Axisymmetric Tokamak Equilibrium in Toroidal Coordinates with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  q₀ = ", equ.q₀)
end


r(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = x[1]
θ(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = x[2]
ϕ(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = x[3]
R(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = equ.R₀ + r(x,equ) * cos(θ(x,equ))
X(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = R(x,equ) * cos(ϕ(x,equ))
Y(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = R(x,equ) * sin(ϕ(x,equ))
Z(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = r(x,equ) * sin(θ(x,equ))

J(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = r(x,equ) * R(x,equ)

A₁(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = zero(eltype(x))
A₂(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = + equ.B₀ * equ.R₀ / cos(θ(x,equ))^2 * ( r(x,equ) * cos(θ(x,equ)) - equ.R₀ * log(R(x,equ) / equ.R₀) )
A₃(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = - equ.B₀ * r(x,equ)^2 / equ.q₀ / 2

x¹(ξ::AbstractVector, equ::AxisymmetricTokamakToroidal) = X(ξ,equ)
x²(ξ::AbstractVector, equ::AxisymmetricTokamakToroidal) = Y(ξ,equ)
x³(ξ::AbstractVector, equ::AxisymmetricTokamakToroidal) = Z(ξ,equ)

ξ¹(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = sqrt((sqrt(x[1]^2 + x[2]^2)-equ.R₀)^2 + x[3]^2)
ξ²(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = atan(x[3], sqrt(x[1]^2 + x[2]^2)-equ.R₀)
ξ³(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = atan(x[2], x[1])

g₁₁(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = one(eltype(x))
g₂₂(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = r(x, equ)^2
g₃₃(x::AbstractVector, equ::AxisymmetricTokamakToroidal) = R(x, equ)^2


get_functions(::AxisymmetricTokamakToroidal) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ)


function periodicity(x::AbstractVector, equ::AxisymmetricTokamakToroidal)
    p = zero(x)
    p[2] = 2π
    p[3] = 2π
    return p
end

macro axisymmetric_tokamak_equilibrium_toroidal(R₀, B₀, q₀)
    generate_equilibrium_code(AxisymmetricTokamakToroidal(R₀, B₀, q₀); output=false)
end
