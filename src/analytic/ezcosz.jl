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
struct EzCosZ{T <: Number} <: CartesianPerturbation
    name::String
    E₀::T
    EzCosZ{T}(E₀::T) where T <: Number = new("EzCosZ", E₀)
end

EzCosZ(E₀::T) where T <: Number = EzCosZ{T}(E₀)


function Base.show(io::IO, equ::EzCosZ)
    print(io, "Simple perturbation in electric field")
end


function φ(x::AbstractArray{T,1}, equ::EzCosZ) where {T <: Number}
   equ.E₀ / (2π) * sin(2π * Z(x,equ))
end


macro ezcosz_perturbation(E₀=1.)
    generate_equilibrium_code(EzCosZ(E₀); output=false)
end
