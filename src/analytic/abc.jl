@doc raw"""
ABCEquilibrium equilibrium in (x,y,z) coordinates with covariant components of the vector
potential given by
```math
A (x,y,z) = \big( a \, \sin(z) + c \, \cos(y) , \, b \, \sin(x) + a \, \cos(z) , \, c \, \sin(y) + b \, \cos(x) \big)^T
```
resulting in the magnetic field ``B(x,y,z) = A(x,y,z)``.

Parameters: `a`, `b`, `c`
"""
module ABC

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..AnalyticCartesianField: X, Y, Z

    export ABCEquilibrium

    struct ABCEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        a::T
        b::T
        c::T

        ABCEquilibrium{T}(a::T, b::T, c::T) where T <: Number = new("ABCEquilibrium", a, b, c)
    end

    ABCEquilibrium(a::T, b::T, c::T) where T <: Number = ABCEquilibrium{T}(a, b, c)


    function Base.show(io::IO, equ::ABCEquilibrium)
        print(io, "ABC Equilibrium with\n")
        print(io, "  A = ", equ.a, "\n")
        print(io, "  B = ", equ.b, "\n")
        print(io, "  C = ", equ.c)
    end


    ElectromagneticFields.A₁(x::AbstractVector, equ::ABCEquilibrium) = equ.a * sin(x[3]) + equ.c * cos(x[2])
    ElectromagneticFields.A₂(x::AbstractVector, equ::ABCEquilibrium) = equ.b * sin(x[1]) + equ.a * cos(x[3])
    ElectromagneticFields.A₃(x::AbstractVector, equ::ABCEquilibrium) = equ.c * sin(x[2]) + equ.b * cos(x[1])

    ElectromagneticFields.get_functions(::ABCEquilibrium) = (X=X, Y=Y, Z=Z)

    macro abc_equilibrium(a, b, c)
        generate_equilibrium_code(ABCEquilibrium(a, b, c); output=false)
    end

    function init(a, b, c; perturbation=ZeroPerturbation())
        equilibrium = ABCEquilibrium(a, b, c)
        load_equilibrium(equilibrium, perturbation; target_module=ABC)
        return equilibrium
    end

end
