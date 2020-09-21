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
module EzCosZ

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..AnalyticCartesianField: X, Y, Z

    export EzCosZPerturbation

    const DEFAULT_E₀ = 1.0

    struct EzCosZPerturbation{T <: Number} <: CartesianPerturbation
        name::String
        E₀::T
        EzCosZPerturbation{T}(E₀::T) where T <: Number = new("EzCosZ", E₀)
    end

    EzCosZPerturbation(E₀::T=DEFAULT_E₀) where T <: Number = EzCosZPerturbation{T}(E₀)


    function Base.show(io::IO, equ::EzCosZPerturbation)
        print(io, "Simple perturbation in electric field")
    end


    ElectromagneticFields.φ(x::AbstractVector, equ::EzCosZPerturbation) = equ.E₀ / (2π) * sin(2π * Z(x,equ))

    ElectromagneticFields.get_functions(::EzCosZPerturbation) = (X=X, Y=Y, Z=Z)


    macro ezcosz_perturbation(E₀=1.)
        generate_equilibrium_code(EzCosZPerturbation(E₀); output=false)
    end

    function init(E₀=DEFAULT_E₀)
        EzCosZPerturbation(E₀)
    end

end
