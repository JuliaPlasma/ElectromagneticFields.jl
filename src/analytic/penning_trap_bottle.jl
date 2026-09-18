
const DEFAULT_PENNING_BOTTLE_B₀ = 100.0
const DEFAULT_PENNING_BOTTLE_Bₚ = 200.0
const DEFAULT_PENNING_BOTTLE_E₀ = 10.0

@doc raw"""
Penning trap with magnetic bottle in (x,y,z) coordinates.
Based on Yanyan Shi, Yajuan Sun, Yulei Wang, Jian Liu, Study of adaptive symplectic methods for
      simulating charged particle dynamics, Journal of Computational Dynamics 6, 429-448, 2019.

The covariant components of the vector potential are given by
```math
A (x,y,z) = B_0 / 2 \, ( -y , x, 0)^T - B_1 \, (xz, yz, (x^2 + y^2)/2 - z^2)^T,
```
resulting in the magnetic field with covariant components
```math
B (x,y,z) = B_0 \, ( 0, 0, 1)^T - B_1 \, ( yz^2 - y^3 / 6, x^3 / 6, xyz )^T ,
```
and the electrostatic potential given by
```math
\varphi (x,y,z) = - E_0 \, ( x^2 / 2 + y^2 / 2 - z^2) ,
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
struct PenningTrapBottleEquilibrium{T <: Number} <: CartesianEquilibrium
    name::String
    B₀::T
    Bₚ::T
    E₀::T

    function PenningTrapBottleEquilibrium{T}(B₀::T, Bₚ::T, E₀::T) where {T <: Number}
        new("PenningTrapBottleEquilibrium", B₀, Bₚ, E₀)
    end
end

function PenningTrapBottleEquilibrium(
        B₀::T = DEFAULT_PENNING_BOTTLE_B₀,
        Bₚ::T = DEFAULT_PENNING_BOTTLE_Bₚ,
        E₀::T = DEFAULT_PENNING_BOTTLE_E₀) where {T <: Number}
    PenningTrapBottleEquilibrium{T}(B₀, Bₚ, E₀)
end

function Base.show(io::IO, equ::PenningTrapBottleEquilibrium)
    print(io, "Penning trap with magnetic bottle in (x,y,z) coordinates with\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  Bₚ = ", equ.Bₚ, "\n")
    print(io, "  E₀ = ", equ.E₀)
end

function A₁(x::AbstractVector, equ::PenningTrapBottleEquilibrium)
    -equ.B₀ / 2 * Y(x, equ) - equ.Bₚ * (Y(x, equ) * Z(x, equ)^2 - Y(x, equ)^3 / 6)
end
function A₂(x::AbstractVector, equ::PenningTrapBottleEquilibrium)
    +equ.B₀ / 2 * X(x, equ) - equ.Bₚ * X(x, equ)^3 / 6
end
function A₃(x::AbstractVector, equ::PenningTrapBottleEquilibrium)
    -equ.Bₚ * X(x, equ) * Y(x, equ) * Z(x, equ)
end

function φ(x::AbstractVector, equ::PenningTrapBottleEquilibrium)
    -equ.E₀ * (X(x, equ)^2 / 2 + Y(x, equ)^2 / 2 - Z(x, equ)^2)
end

get_functions(::PenningTrapBottleEquilibrium) = (X = X, Y = Y, Z = Z)
