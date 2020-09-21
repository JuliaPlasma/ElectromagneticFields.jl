@doc raw"""
θ-pinch equilibrium in (x,y,z) coordinates with covariant components of
the vector potential given by
```math
A (x,y,z) = \frac{B_0}{2} \big( - y , \, x , \, 0 \big)^T
```
resulting in the magnetic field with covariant components
```math
B (x,y,z) = \big( 0 , \, 0 , \, B_0 \big)^T .
```

Parameters:
    B₀: B-field at magnetic axis
"""
module ThetaPinch

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..AnalyticCartesianField: X, Y, Z

    export ThetaPinchEquilibrium

    const DEFAULT_B₀ = 1.0

    struct ThetaPinchEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        B₀::T

        function ThetaPinchEquilibrium{T}(B₀::T) where T <: Number
            new("ThetaPinchEquilibrium", B₀)
        end
    end

    ThetaPinchEquilibrium(B₀::T=DEFAULT_B₀) where T <: Number = ThetaPinchEquilibrium{T}(B₀)


    function Base.show(io::IO, equ::ThetaPinchEquilibrium)
        print(io, "θ-Pinch Equilibrium in (x,y,z) Coordinates with\n")
        print(io, "  B₀ = ", equ.B₀)
    end


    r²(x::AbstractVector, equ::ThetaPinchEquilibrium) = X(x,equ)^2 + Y(x,equ)^2
    r(x::AbstractVector, equ::ThetaPinchEquilibrium) = sqrt(r²(x,equ))
    R(x::AbstractVector, equ::ThetaPinchEquilibrium) = r(x,equ)
    θ(x::AbstractVector, equ::ThetaPinchEquilibrium) = atan(Y(x,equ), X(x,equ))
    ϕ(x::AbstractVector, equ::ThetaPinchEquilibrium) = θ(x,equ)

    ElectromagneticFields.A₁(x::AbstractVector, equ::ThetaPinchEquilibrium) = - equ.B₀ * Y(x,equ) / 2
    ElectromagneticFields.A₂(x::AbstractVector, equ::ThetaPinchEquilibrium) = + equ.B₀ * X(x,equ) / 2
    ElectromagneticFields.A₃(x::AbstractVector, equ::ThetaPinchEquilibrium) = zero(eltype(x))

    ElectromagneticFields.get_functions(::ThetaPinchEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)


    macro zpinch_equilibrium(B₀=DEFAULT_B₀)
        generate_equilibrium_code(ThetaPinchEquilibrium(B₀); output=false)
    end

    function init(B₀=DEFAULT_B₀; perturbation=ZeroPerturbation())
        equilibrium = ThetaPinchEquilibrium(B₀)
        load_equilibrium(equilibrium, perturbation; target_module=ThetaPinch)
        return equilibrium
    end

end
