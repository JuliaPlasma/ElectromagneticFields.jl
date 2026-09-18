
const DEFAULT_EZCOSZ_E₀ = 1.0

@doc raw"""
Simple perturbation in electric field in (x,y,z) coordinates:
```math
E(x,y,z) = E_0 \, \begin{pmatrix}
0 \\
0 \\
\cos (2\pi z) \\
\end{pmatrix}
```

Parameters: `E₀`
"""
struct EzCosZPerturbation{T <: Number} <: CartesianPerturbation
    name::String
    E₀::T
    EzCosZPerturbation{T}(E₀::T) where {T <: Number} = new("EzCosZ", E₀)
end

function EzCosZPerturbation(E₀::T = DEFAULT_EZCOSZ_E₀) where {T <: Number}
    EzCosZPerturbation{T}(E₀)
end

function Base.show(io::IO, equ::EzCosZPerturbation)
    print(io, "Simple perturbation in electric field")
end

φ(x::AbstractVector, equ::EzCosZPerturbation) = equ.E₀ / (2π) * sin(2π * Z(x, equ))

get_functions(::EzCosZPerturbation) = (X = X, Y = Y, Z = Z)
