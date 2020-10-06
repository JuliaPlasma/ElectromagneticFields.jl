"""
Axisymmetric Solov'ev equilibra in (R/R₀,Z/R₀,phi) coordinates.
Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `ϵ`:  inverse aspect ratio
 * `κ`:  elongation
 * `δ`:  triangularity
 * `α`:  free constant, determined to match a given beta value
"""
module Solovev

    using RecipesBase

    using SymPy: N, Sym, diff, expand, solve, subs

    import ..ElectromagneticFields
    import ..ElectromagneticFields: ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..SolovevAbstract: AbstractSolovevEquilibrium, X, Y, Z, R, r, θ, ϕ, r²

    export  SolovevEquilibrium

    include("solovev_psi.jl")

    struct SolovevEquilibrium{T <: Number} <: AbstractSolovevEquilibrium
        name::String
        R₀::T
        B₀::T
        ϵ::T
        κ::T
        δ::T
        α::T
        c::Vector{T}

        function SolovevEquilibrium{T}(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, α::T, c::Vector{T}) where T <: Number
            new("Solovev Equilibrium", R₀, B₀, ϵ, κ, δ, α, c)
        end
    end

    function SolovevEquilibrium(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, α::T) where T <: Number

        n = 7

        x = [Sym("x" * string(i)) for i in 1:3]
        c = [Sym("c" * string(i)) for i in 1:n]

        ψ = ( ψ₀(x,α) + c[1] * ψ₁(x)
                      + c[2] * ψ₂(x)
                      + c[3] * ψ₃(x)
                      + c[4] * ψ₄(x)
                      + c[5] * ψ₅(x)
                      + c[6] * ψ₆(x)
                      + c[7] * ψ₇(x) )

        eqs = [
            expand(subs(subs(ψ, x[1], 1+ϵ), x[2], 0)),
            expand(subs(subs(ψ, x[1], 1-ϵ), x[2], 0)),
            expand(subs(subs(ψ, x[1], 1-δ*ϵ), x[2], κ*ϵ)),
            expand(subs(subs(diff(ψ, x[1]), x[1], 1-δ*ϵ), x[2], κ*ϵ)),
            expand(subs(subs(diff(ψ, x[2], 2), x[1], 1+ϵ), x[2], 0) - (1 + asin(δ))^2 / (ϵ * κ^2) * subs(subs(diff(ψ, x[1]), x[1], 1+ϵ), x[2], 0)),
            expand(subs(subs(diff(ψ, x[2], 2), x[1], 1-ϵ), x[2], 0) + (1 - asin(δ))^2 / (ϵ * κ^2) * subs(subs(diff(ψ, x[1]), x[1], 1-ϵ), x[2], 0)),
            expand(subs(subs(diff(ψ, x[1], 2), x[1], 1-δ*ϵ), x[2], κ*ϵ) - κ / (ϵ * (1 - δ^2)) * subs(subs(diff(ψ, x[2]), x[1], 1-δ*ϵ), x[2], κ*ϵ))
        ]

        csym = solve(eqs, c)
        cnum = [N(csym[c[i]]) for i in 1:n]

        SolovevEquilibrium{T}(R₀, B₀, ϵ, κ, δ, α, cnum)
    end


    SolovevEquilibriumITER() = SolovevEquilibrium(6.2, 5.3, 0.32, 1.7, 0.33, -0.155)
    SolovevEquilibriumNSTX() = SolovevEquilibrium(0.85, 0.3, 0.78, 2.0, 0.35, 1.0)
    SolovevEquilibriumFRC()  = SolovevEquilibrium(0.0, 0.0, 0.99, 10., 0.7, 0.0)
    # SolovevEquilibriumFRC2() = SolovevEquilibrium(0.0, 0.0, 1.00, 10., 1.0, 0.0)


    function Base.show(io::IO, equ::SolovevEquilibrium)
        print(io, "SolovevEquilibrium Equilibrium with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  ϵ  = ", equ.ϵ,  "\n")
        print(io, "  κ  = ", equ.κ,  "\n")
        print(io, "  δ  = ", equ.δ,  "\n")
        print(io, "  α  = ", equ.α)
    end


    function ElectromagneticFields.A₃(x::AbstractArray{T,1}, equ::SolovevEquilibrium) where {T <: Number}
        ( ψ₀(x, equ.α) + equ.c[1] * ψ₁(x)
                       + equ.c[2] * ψ₂(x)
                       + equ.c[3] * ψ₃(x)
                       + equ.c[4] * ψ₄(x)
                       + equ.c[5] * ψ₅(x)
                       + equ.c[6] * ψ₆(x)
                       + equ.c[7] * ψ₇(x) )
    end


    macro solovev_equilibrium(R₀, B₀, ϵ, κ, δ, α)
        generate_equilibrium_code(SolovevEquilibrium(R₀, B₀, ϵ, κ, δ, α); output=false)
    end

    function init(R₀, B₀, ϵ, κ, δ, α; perturbation=ZeroPerturbation(), target_module=Solovev)
        equilibrium = SolovevEquilibrium(R₀, B₀, ϵ, κ, δ, α)
        load_equilibrium(equilibrium, perturbation; target_module=target_module)
        return equilibrium
    end

    function ITER(; perturbation=ZeroPerturbation(), target_module=Solovev)
        equilibrium = SolovevEquilibriumITER()
        load_equilibrium(equilibrium, perturbation; target_module=target_module)
        return equilibrium
    end

    function NSTX(; perturbation=ZeroPerturbation(), target_module=Solovev)
        equilibrium = SolovevEquilibriumNSTX()
        load_equilibrium(equilibrium, perturbation; target_module=target_module)
        return equilibrium
    end

    function FRC(; perturbation=ZeroPerturbation(), target_module=Solovev)
        equilibrium = SolovevEquilibriumFRC()
        load_equilibrium(equilibrium, perturbation; target_module=target_module)
        return equilibrium
    end


    @recipe function f(equ::SolovevEquilibrium;
                       nx = 100, ny = 120, nτ = 200, levels = 50, size = (300,400), aspect_ratio = :equal,
                       xlims = ( 0.50,  1.50),
                       ylims = (-0.75, +0.75))

        xgrid = LinRange(xlims[1], xlims[2], nx)
        zgrid = LinRange(ylims[1], ylims[2], ny)
        pot   = [A₃(0, xgrid[i], zgrid[j], 0.0) / xgrid[i] for i in eachindex(xgrid), j in eachindex(zgrid)]

        τ = LinRange(0, 2π, nτ)
        boundary_X = 1 .+ equ.ϵ .* cos.(τ .+ asin(equ.δ) .* sin.(τ) )
        boundary_Y = equ.ϵ .* equ.κ .* sin.(τ)

        aspect_ratio := aspect_ratio
        size   := size
        xlims  := xlims
        ylims  := ylims
        levels := levels
        legend := :none

        @series begin
            seriestype := :contour
            (xgrid, zgrid, pot')
        end

        @series begin
            seriestype  := :path
            seriescolor := :red
            linewidth := 3
            (boundary_X, boundary_Y)
        end
    end
    
end
