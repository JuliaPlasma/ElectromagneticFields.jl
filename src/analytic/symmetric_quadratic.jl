@doc raw"""
Symmetric quadratic mangetic field in (x,y,z) coordinates with covariant components of
the vector potential given by
```math
A (x,y,z) = \frac{B_0}{4} \big( - y \, (2 + x^2 + y^2) , \, x \, (2 + x^2 + y^2) , \, 0 \big)^T
```
resulting in the magnetic field with covariant components
```math
B(x,y,z) = B_0 \, \begin{pmatrix}
0 \\
0 \\
1 + x^2 + y^2 \\
\end{pmatrix}
```

Parameters: `B₀`
"""
module SymmetricQuadratic

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..AnalyticCartesianField: X, Y, Z

    export SymmetricQuadraticEquilibrium

    const DEFAULT_B₀ = 1.0

    struct SymmetricQuadraticEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        B₀::T
        SymmetricQuadraticEquilibrium{T}(B₀::T) where T <: Number = new("Symmetric Quadratic Magnetic Field", B₀)
    end

    SymmetricQuadraticEquilibrium(B₀::T=DEFAULT_B₀) where T <: Number = SymmetricQuadraticEquilibrium{T}(B₀)


    function Base.show(io::IO, equ::SymmetricQuadraticEquilibrium)
        print(io, equ.name)
    end


    r²(x::AbstractVector, equ::SymmetricQuadraticEquilibrium) = X(x,equ)^2 + Y(x,equ)^2
    r(x::AbstractVector, equ::SymmetricQuadraticEquilibrium) = sqrt(r²(x,equ))
    R(x::AbstractVector, equ::SymmetricQuadraticEquilibrium) = r(x,equ)
    θ(x::AbstractVector, equ::SymmetricQuadraticEquilibrium) = atan(Y(x,equ), X(x,equ))
    ϕ(x::AbstractVector, equ::SymmetricQuadraticEquilibrium) = θ(x,equ)

    ElectromagneticFields.A₁(x::AbstractVector, equ::SymmetricQuadraticEquilibrium) = - equ.B₀ * x[2] * (2 + x[1]^2 + x[2]^2) / 4
    ElectromagneticFields.A₂(x::AbstractVector, equ::SymmetricQuadraticEquilibrium) = + equ.B₀ * x[1] * (2 + x[1]^2 + x[2]^2) / 4
    ElectromagneticFields.A₃(x::AbstractVector, equ::SymmetricQuadraticEquilibrium) = zero(eltype(x))

    ElectromagneticFields.get_functions(::SymmetricQuadraticEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)

    macro symmetric_quadratic_equilibrium(B₀=DEFAULT_B₀)
        generate_equilibrium_code(SymmetricQuadraticEquilibrium(B₀); output=false)
    end

    function init(B₀=DEFAULT_B₀; perturbation=ZeroPerturbation())
        equilibrium = SymmetricQuadraticEquilibrium(B₀)
        load_equilibrium(equilibrium, perturbation; target_module=SymmetricQuadratic)
        return equilibrium
    end

end
