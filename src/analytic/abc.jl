@doc raw"""
Arnold-Beltrami-Childress (ABC) field in (x,y,z) coordinates with covariant components of the vector
potential given by
```math
A (x,y,z) = \big( a \, \sin(z) + c \, \cos(y) , \, b \, \sin(x) + a \, \cos(z) , \, c \, \sin(y) + b \, \cos(x) \big)^T
```
resulting in the magnetic field ``B(x,y,z) = A(x,y,z)``.

Parameters: `a`, `b`, `c`
"""
module ABC

    using RecipesBase

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, code
    import ..AnalyticCartesianField: X, Y, Z

    export ABCEquilibrium

    const DEFAULT_A = 1
    const DEFAULT_B = 1
    const DEFAULT_C = 1

    struct ABCEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        a₀::T
        b₀::T
        c₀::T

        ABCEquilibrium{T}(a::T, b::T, c::T) where T <: Number = new("ABCEquilibrium", a, b, c)
    end

    ABCEquilibrium(a::T=DEFAULT_A, b::T=DEFAULT_B, c::T=DEFAULT_C) where T <: Number = ABCEquilibrium{T}(a, b, c)

    function init(a=DEFAULT_A, b=DEFAULT_B, c=DEFAULT_C)
        ABCEquilibrium(a, b, c)
    end

    macro code(a=DEFAULT_A, b=DEFAULT_B, c=DEFAULT_C)
        code(init(a, b, c); escape=true)
    end

    function Base.show(io::IO, equ::ABCEquilibrium)
        print(io, "ABC Equilibrium with\n")
        print(io, "  A = ", equ.a₀, "\n")
        print(io, "  B = ", equ.b₀, "\n")
        print(io, "  C = ", equ.c₀)
    end


    ElectromagneticFields.A₁(x::AbstractVector, equ::ABCEquilibrium) = equ.a₀ * sin(x[3]) + equ.c₀ * cos(x[2])
    ElectromagneticFields.A₂(x::AbstractVector, equ::ABCEquilibrium) = equ.b₀ * sin(x[1]) + equ.a₀ * cos(x[3])
    ElectromagneticFields.A₃(x::AbstractVector, equ::ABCEquilibrium) = equ.c₀ * sin(x[2]) + equ.b₀ * cos(x[1])

    ElectromagneticFields.get_functions(::ABCEquilibrium) = (X=X, Y=Y, Z=Z)


    @recipe function f(equ::ABCEquilibrium; nx=99, ni=div(nx,2)+1, nl=12, size=(400,1200))
        @eval $(code(equ)) # this causes world age problems

        xmin = 0
        xmax = 2π
        grid = LinRange(xmin, xmax, nx)
        Bfield = [B(0, grid[i], grid[j], grid[k]) for i in eachindex(grid), j in eachindex(grid), k in eachindex(grid)]

        seriestype   := :contour
        aspect_ratio := :equal
        layout := (3, 1)
        size   := size
        xlims  := (xmin,xmax)
        ylims  := (xmin,xmax)
        levels := nl
        legend := :none

        @series begin
            subplot := 1
            xguide := "x"
            yguide := "y"
            title  := "|B(x,y,π)|"
            (grid, grid, Bfield[:,:,ni])
        end

        @series begin
            subplot := 2
            xguide := "x"
            yguide := "z"
            title  := "|B(x,π,z)|"
            (grid, grid, Bfield[:,ni,:])
        end

        @series begin
            subplot := 3
            xguide := "y"
            yguide := "z"
            title  := "|B(π,y,z)|"
            (grid, grid, Bfield[ni,:,:])
        end
    end

end
