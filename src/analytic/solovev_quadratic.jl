"""
Quadratic Solov'ev equilibrium in (R,Z,phi) coordinates.
Based on McCarthy, Physics of Plasmas 6, 3554, 1999.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `a`, b`: free constants
"""
module SolovevQuadratic

    import ..ElectromagneticFields
    import ..ElectromagneticFields: ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..SolovevAbstract: AbstractSolovevEquilibrium, X, Y, Z, R, r, θ, ϕ, r²

    export  SolovevQuadraticEquilibrium

    struct SolovevQuadraticEquilibrium{T <: Number} <: AbstractSolovevEquilibrium
        name::String
        R₀::T
        B₀::T
        a::T
        b::T

        function SolovevQuadraticEquilibrium{T}(R₀::T, B₀::T, a::T, b::T) where T <: Number
            new("QuadraticSolovevEquilibrium", R₀, B₀, a, b)
        end
    end

    function SolovevQuadraticEquilibrium(R₀::T, B₀::T, a::T, b::T) where T <: Number
        SolovevQuadraticEquilibrium{T}(R₀, B₀, a, b)
    end


    function Base.show(io::IO, equ::SolovevQuadraticEquilibrium)
        print(io, "Quadratic Solovev Equilibrium with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  a  = ", equ.a, "\n")
        print(io, "  b  = ", equ.b)
    end


    ElectromagneticFields.A₁(x::AbstractVector, equ::SolovevQuadraticEquilibrium) = zero(eltype(x))
    ElectromagneticFields.A₂(x::AbstractVector, equ::SolovevQuadraticEquilibrium) = zero(eltype(x))
    ElectromagneticFields.A₃(x::AbstractVector, equ::SolovevQuadraticEquilibrium) = - equ.a * R(x,equ)^4 / 8 - equ.b * Z(x,equ)^2 / 2


    macro solovev_equilibrium_quadratic(R₀, B₀, a, b)
        generate_equilibrium_code(SolovevQuadraticEquilibrium(R₀, B₀, a, b); output=false)
    end

    function init(R₀, B₀, a, b; perturbation=ZeroPerturbation())
        equilibrium = SolovevQuadraticEquilibrium(R₀, B₀, a, b)
        load_equilibrium(equilibrium, perturbation; target_module=SolovevQuadratic)
        return equilibrium
    end

end
