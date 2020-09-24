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
module AxisymmetricTokamakCircular

    import ..ElectromagneticFields
    import ..ElectromagneticFields: AnalyticEquilibrium, ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code

    export  AxisymmetricTokamakCircularEquilibrium

    const DEFAULT_R₀ = 1.0
    const DEFAULT_B₀ = 1.0
    const DEFAULT_q₀ = 2.0

    struct AxisymmetricTokamakCircularEquilibrium{T <: Number} <: AnalyticEquilibrium
        name::String
        R₀::T
        B₀::T
        q₀::T

        function AxisymmetricTokamakCircularEquilibrium{T}(R₀::T, B₀::T, q₀::T) where T <: Number
            new("AxisymmetricTokamakEquilibriumToroidal", R₀, B₀, q₀)
        end
    end

    AxisymmetricTokamakCircularEquilibrium(R₀::T=DEFAULT_R₀, B₀::T=DEFAULT_B₀, q₀::T=DEFAULT_q₀) where T <: Number = AxisymmetricTokamakCircularEquilibrium{T}(R₀, B₀, q₀)

    function Base.show(io::IO, equ::AxisymmetricTokamakCircularEquilibrium)
        print(io, "Axisymmetric Tokamak Equilibrium in Toroidal Coordinates with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  q₀ = ", equ.q₀)
    end


    r(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = x[1]
    θ(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = x[2]
    ϕ(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = x[3]
    R(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = equ.R₀ + r(x,equ) * cos(θ(x,equ))
    X(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = R(x,equ) * cos(ϕ(x,equ))
    Y(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = R(x,equ) * sin(ϕ(x,equ))
    Z(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = r(x,equ) * sin(θ(x,equ))

    ElectromagneticFields.J(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = r(x,equ) * R(x,equ)
    
    ElectromagneticFields.A₁(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = + equ.B₀ * equ.R₀ * ( Z(x,equ) / R(x,equ) * cos(θ(x,equ)) - log(R(x,equ) / equ.R₀) * sin(θ(x,equ)) ) / 2
    ElectromagneticFields.A₂(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = - equ.B₀ * equ.R₀ * ( Z(x,equ) / R(x,equ) * sin(θ(x,equ)) + log(R(x,equ) / equ.R₀) * cos(θ(x,equ)) ) * r(x,equ) / 2
    ElectromagneticFields.A₃(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = + equ.B₀ * r(x,equ)^2 / equ.q₀ / 2

    ElectromagneticFields.x¹(ξ::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = X(ξ,equ)
    ElectromagneticFields.x²(ξ::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = Y(ξ,equ)
    ElectromagneticFields.x³(ξ::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = Z(ξ,equ)

    ElectromagneticFields.ξ¹(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = sqrt((sqrt(x[1]^2 + x[2]^2)-equ.R₀)^2 + x[3]^2)
    ElectromagneticFields.ξ²(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = atan(x[3], sqrt(x[1]^2 + x[2]^2)-equ.R₀)
    ElectromagneticFields.ξ³(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = atan(x[2], x[1])

    ElectromagneticFields.g₁₁(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = one(eltype(x))
    ElectromagneticFields.g₂₂(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = r(x, equ)^2
    ElectromagneticFields.g₃₃(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium) = R(x, equ)^2

    ElectromagneticFields.get_functions(::AxisymmetricTokamakCircularEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ)

    function ElectromagneticFields.periodicity(x::AbstractVector, equ::AxisymmetricTokamakCircularEquilibrium)
        p = zero(x)
        p[2] = 2π
        p[3] = 2π
        return p
    end

    macro axisymmetric_tokamak_equilibrium_toroidal(R₀, B₀, q₀)
        generate_equilibrium_code(AxisymmetricTokamakCircularEquilibrium(R₀, B₀, q₀); output=false)
    end

    function init(R₀=DEFAULT_R₀, B₀=DEFAULT_B₀, q₀=DEFAULT_q₀; perturbation=ZeroPerturbation())
        equilibrium = AxisymmetricTokamakCircularEquilibrium(R₀, B₀, q₀)
        load_equilibrium(equilibrium, perturbation; target_module=AxisymmetricTokamakCircular)
        return equilibrium
    end

end
