@doc raw"""
Singular magnetic field in (x,y,z) coordinates with covariant components of
the vector potential given by
```math
A (x,y,z) = \frac{B_0}{\sqrt{(x^2 + y^2)}^3} \big( y , \, - x , \, 0 \big)^T
```
resulting in the magnetic field with covariant components
```math
B(x,y,z) = B_0 \, \begin{pmatrix}
0 \\
0 \\
(x^2 + y^2)^{-3/2} \\
\end{pmatrix}
```

Parameters: `B₀`
"""
module Singular

    using RecipesBase
    using LaTeXStrings

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..AnalyticCartesianField: X, Y, Z

    export SingularEquilibrium

    const DEFAULT_B₀ = 1.0

    struct SingularEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        B₀::T
        SingularEquilibrium{T}(B₀::T) where T <: Number = new("Singular Magnetic Field", B₀)
    end

    SingularEquilibrium(B₀::T=DEFAULT_B₀) where T <: Number = SingularEquilibrium{T}(B₀)


    function Base.show(io::IO, equ::SingularEquilibrium)
        print(io, equ.name)
    end


    r²(x::AbstractVector, equ::SingularEquilibrium) = X(x,equ)^2 + Y(x,equ)^2
    r(x::AbstractVector, equ::SingularEquilibrium) = sqrt(r²(x,equ))
    R(x::AbstractVector, equ::SingularEquilibrium) = r(x,equ)
    θ(x::AbstractVector, equ::SingularEquilibrium) = atan(Y(x,equ), X(x,equ))
    ϕ(x::AbstractVector, equ::SingularEquilibrium) = θ(x,equ)

    ElectromagneticFields.A₁(x::AbstractVector, equ::SingularEquilibrium) = - equ.B₀ * x[2] * (2 + x[1]^2 + x[2]^2) / 4
    ElectromagneticFields.A₂(x::AbstractVector, equ::SingularEquilibrium) = + equ.B₀ * x[1] * (2 + x[1]^2 + x[2]^2) / 4
    ElectromagneticFields.A₃(x::AbstractVector, equ::SingularEquilibrium) = zero(eltype(x))

    ElectromagneticFields.get_functions(::SingularEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)

    macro symmetric_quadratic_equilibrium(B₀=DEFAULT_B₀)
        generate_equilibrium_code(SingularEquilibrium(B₀); output=false)
    end

    function init(B₀=DEFAULT_B₀; perturbation=ZeroPerturbation())
        equilibrium = SingularEquilibrium(B₀)
        load_equilibrium(equilibrium, perturbation; target_module=Singular)
        return equilibrium
    end


    @recipe function f(equ::SingularEquilibrium;
                       nx = 100, ny = 100, levels = 20, size = (400,1200),
                       xlims = (-1., +1.),
                       ylims = (-1., +1.))

        xgrid = LinRange(xlims[1], xlims[2], nx)
        ygrid = LinRange(ylims[1], ylims[2], ny)
        pot1  = [A₁(0, xgrid[i], ygrid[j], 0.0) for i in eachindex(xgrid), j in eachindex(ygrid)]
        pot2  = [A₂(0, xgrid[i], ygrid[j], 0.0) for i in eachindex(xgrid), j in eachindex(ygrid)]
        Bfield = [B₃(0, xgrid[i], ygrid[j], 0.0) for i in eachindex(xgrid), j in eachindex(ygrid)]

        seriestype   := :contour
        aspect_ratio := :equal
        layout := (3,1)
        size   := size
        xlims  := xlims
        ylims  := ylims
        levels := levels
        legend := :none

        @series begin
            subplot := 1
            title  := L"A_x (x,y)"
            xguide := L"x"
            yguide := L"y"
            (xgrid, ygrid, pot1)
        end

        @series begin
            subplot := 2
            title  := L"A_y (x,y)"
            xguide := L"x"
            yguide := L"y"
            (xgrid, ygrid, pot2)
        end

        @series begin
            subplot := 3
            title  := L"B_z (x,y)"
            xguide := L"x"
            yguide := L"y"
            (xgrid, ygrid, Bfield)
        end
    end

end
