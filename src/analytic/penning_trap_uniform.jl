@doc raw"""
Penning trap with uniform magnetic field in (x,y,z) coordinates.
Based on Yanyan Shi, Yajuan Sun, Yulei Wang, Jian Liu, Study of adaptive symplectic methods for
      simulating charged particle dynamics, Journal of Computational Dynamics 6, 429-448, 2019.

The covariant components of the vector potential are given by
```math
A (x,y,z) = B_0 \, ( 0, x, 0)^T ,
```
resulting in the magnetic field with covariant components
```math
B (x,y,z) = B_0 \, ( 0, 0, 1)^T ,
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
* `E₀`: E-field strength
"""
module PenningTrapUniform

    using RecipesBase

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, code
    import ..AnalyticCartesianField: X, Y, Z

    export  PenningTrapUniformEquilibrium

    const DEFAULT_B₀ = 100.0
    const DEFAULT_E₀ = 10.0

    struct PenningTrapUniformEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        B₀::T
        E₀::T

        function PenningTrapUniformEquilibrium{T}(B₀::T, E₀::T) where T <: Number
            new("PenningTrapUniformEquilibrium", B₀, E₀)
        end
    end

    PenningTrapUniformEquilibrium(B₀::T=DEFAULT_B₀, E₀::T=DEFAULT_E₀) where T <: Number = PenningTrapUniformEquilibrium{T}(B₀, E₀)

    function init(B₀=DEFAULT_B₀, E₀=DEFAULT_E₀)
        PenningTrapUniformEquilibrium(B₀, E₀)
    end
    
    macro code(B₀=DEFAULT_B₀, E₀=DEFAULT_E₀)
        code(init(B₀, E₀); escape=true)
    end

    function Base.show(io::IO, equ::PenningTrapUniformEquilibrium)
        print(io, "Penning trap with uniform magnetic field in (x,y,z) coordinates with\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  E₀ = ", equ.E₀)
    end


    ElectromagneticFields.A₁(x::AbstractVector, equ::PenningTrapUniformEquilibrium) = zero(eltype(x))
    ElectromagneticFields.A₂(x::AbstractVector, equ::PenningTrapUniformEquilibrium) = equ.B₀ * X(x,equ)
    ElectromagneticFields.A₃(x::AbstractVector, equ::PenningTrapUniformEquilibrium) = zero(eltype(x))
    ElectromagneticFields.φ(x::AbstractVector, equ::PenningTrapUniformEquilibrium) = - equ.E₀ * (X(x,equ)^2 / 2 + Y(x,equ)^2 / 2 - Z(x,equ)^2)

    ElectromagneticFields.get_functions(::PenningTrapUniformEquilibrium) = (X=X, Y=Y, Z=Z)

end
