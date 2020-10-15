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

    using LaTeXStrings
    # using Plots
    using RecipesBase

    import ..ElectromagneticFields
    import ..ElectromagneticFields: CartesianEquilibrium, code
    import ..AnalyticCartesianField: X, Y, Z

    export SingularEquilibrium

    const DEFAULT_B₀ = 1.0

    struct SingularEquilibrium{T <: Number} <: CartesianEquilibrium
        name::String
        B₀::T
        SingularEquilibrium{T}(B₀::T) where T <: Number = new("Singular Magnetic Field", B₀)
    end

    SingularEquilibrium(B₀::T=DEFAULT_B₀) where T <: Number = SingularEquilibrium{T}(B₀)

    function init(B₀=DEFAULT_B₀)
        SingularEquilibrium(B₀)
    end

    macro code(B₀=DEFAULT_B₀)
        code(init(B₀); escape=true)
    end


    function Base.show(io::IO, equ::SingularEquilibrium)
        print(io, equ.name)
    end


    r²(x::AbstractVector, equ::SingularEquilibrium) = X(x,equ)^2 + Y(x,equ)^2
    r(x::AbstractVector, equ::SingularEquilibrium) = sqrt(r²(x,equ))
    R(x::AbstractVector, equ::SingularEquilibrium) = r(x,equ)
    θ(x::AbstractVector, equ::SingularEquilibrium) = atan(Y(x,equ), X(x,equ))
    ϕ(x::AbstractVector, equ::SingularEquilibrium) = θ(x,equ)

    ElectromagneticFields.A₁(x::AbstractVector, equ::SingularEquilibrium) = + equ.B₀ * x[2] / sqrt(x[1]^2 + x[2]^2)^3
    ElectromagneticFields.A₂(x::AbstractVector, equ::SingularEquilibrium) = - equ.B₀ * x[1] / sqrt(x[1]^2 + x[2]^2)^3
    ElectromagneticFields.A₃(x::AbstractVector, equ::SingularEquilibrium) = zero(eltype(x))

    ElectromagneticFields.get_functions(::SingularEquilibrium) = (X=X, Y=Y, Z=Z, R=R, r=r, θ=θ, ϕ=ϕ, r²=r²)


    @recipe function f(equ::SingularEquilibrium;
                       nx = 100, ny = 100, levels = 25, size = (400,1200),
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

        logrange(x1, x2, n) = collect(10^y for y in range(log10(x1), log10(x2), length=n))

        function doublelogrange(x1, x2, n)
            lvls = logrange(x1, x2, n)
            vcat(-lvls, +lvls)
        end

        @series begin
            subplot := 1
            title  := L"A_x (x,y)"
            xguide := L"x"
            yguide := L"y"
            levels := doublelogrange(0.1, maximum(pot1), levels)
            # seriescolor := cgrad(:default, levels, scale = :log)
            (xgrid, ygrid, pot1)
        end

        @series begin
            subplot := 2
            title  := L"A_y (x,y)"
            xguide := L"x"
            yguide := L"y"
            levels := doublelogrange(0.1, maximum(pot2), levels)
            # seriescolor := cgrad(:default, levels, scale = :log)
            (xgrid, ygrid, pot2)
        end

        @series begin
            subplot := 3
            title  := L"B_z (x,y)"
            xguide := L"x"
            yguide := L"y"
            levels := logrange(maximum([0.1, minimum(Bfield)]), maximum(Bfield), levels)
            # seriescolor := cgrad(:default, levels, scale = :log)
            (xgrid, ygrid, Bfield)
        end
    end

end
