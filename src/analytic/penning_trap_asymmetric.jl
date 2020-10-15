@doc raw"""
Penning trap with asymmetric magnetic field in (x,y,z) coordinates.
Based on Yanyan Shi, Yajuan Sun, Yulei Wang, Jian Liu, Study of adaptive symplectic methods for
      simulating charged particle dynamics, Journal of Computational Dynamics 6, 429-448, 2019.

The covariant components of the vector potential are given by
```math
A (x,y,z) = B_0 / 2 \, ( -y , x-z/6, y / 6)^T + B_1 / 2 \, (z^2 - y^2, z^2 - x^2, y^2 - x^2)^T,
```
resulting in the magnetic field with covariant components
```math
B (x,y,z) = B_0 \, ( 1/3, 0, 1)^T + B_1 \, ( y-z, x+z, y-x )^T ,
```
and the electrostatic potential given by
```math
\varphi (x,y,z) = E_0 \, ( x^2 / 2 + y^2 / 2 - z^2) ,
```
resulting in the electric field with covariant components
```math
E (x,y,z) = E_0 \, ( x, y, - 2 z)^T .
```

Parameters:
* `B₀`: B-field strength
* `Bₚ`: B-field perturbation strength
* `E₀`: E-field strength
"""
module PenningTrapAsymmetric

    using RecipesBase

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, code
    import ..AnalyticCartesianField: X, Y, Z

    export  PenningTrapAsymmetricEquilibrium

    const DEFAULT_B₀ = 100.0
    const DEFAULT_Bₚ = 50.0
    const DEFAULT_E₀ = 10.0

    struct PenningTrapAsymmetricEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        B₀::T
        Bₚ::T
        E₀::T

        function PenningTrapAsymmetricEquilibrium{T}(B₀::T, Bₚ::T, E₀::T) where T <: Number
            new("PenningTrapAsymmetricEquilibrium", B₀, Bₚ, E₀)
        end
    end

    PenningTrapAsymmetricEquilibrium(B₀::T=DEFAULT_B₀, Bₚ::T=DEFAULT_Bₚ, E₀::T=DEFAULT_E₀) where T <: Number = PenningTrapAsymmetricEquilibrium{T}(B₀, Bₚ, E₀)

    function init(B₀=DEFAULT_B₀, Bₚ=DEFAULT_Bₚ, E₀=DEFAULT_E₀)
        PenningTrapAsymmetricEquilibrium(B₀, E₀, Bₚ)
    end

    macro code(B₀=DEFAULT_B₀, Bₚ=DEFAULT_Bₚ, E₀=DEFAULT_E₀)
        code(init(B₀, Bₚ, E₀); escape=true)
    end

    function Base.show(io::IO, equ::PenningTrapAsymmetricEquilibrium)
        print(io, "Penning trap with asymmetric magnetic field in (x,y,z) coordinates with\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  Bₚ = ", equ.Bₚ, "\n")
        print(io, "  E₀ = ", equ.E₀)
    end


    ElectromagneticFields.A₁(x::AbstractVector, equ::PenningTrapAsymmetricEquilibrium) = - equ.B₀ / 2 *  Y(x,equ)                  + equ.Bₚ / 2* (Z(x,equ)^2 - Y(x,equ)^2)
    ElectromagneticFields.A₂(x::AbstractVector, equ::PenningTrapAsymmetricEquilibrium) = + equ.B₀ / 2 * (X(x,equ) - Z(x,equ) / 6 ) + equ.Bₚ / 2* (Z(x,equ)^2 - X(x,equ)^2)
    ElectromagneticFields.A₃(x::AbstractVector, equ::PenningTrapAsymmetricEquilibrium) = + equ.B₀ / 2 *             Y(x,equ) / 6   + equ.Bₚ / 2* (Y(x,equ)^2 - X(x,equ)^2)
    ElectromagneticFields.φ(x::AbstractVector, equ::PenningTrapAsymmetricEquilibrium) = - equ.E₀ * (X(x,equ)^2 / 2 + Y(x,equ)^2 / 2 - Z(x,equ)^2)

    ElectromagneticFields.get_functions(::PenningTrapAsymmetricEquilibrium) = (X=X, Y=Y, Z=Z)

end
