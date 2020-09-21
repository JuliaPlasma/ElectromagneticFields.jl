"""
Axisymmetric Solov'ev equilibra in (R/R₀,Z/R₀,phi) coordinates.
Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010.

Parameters:
 * `R₀`: position of magnetic axis
 * `B₀`: B-field at magnetic axis
 * `ϵ`:  inverse aspect ratio
 * `κ`:  elongation
 * `δ`:  triangularity
 * `a`:  free constant, determined to match a given beta value
"""
module Solovev

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
        a::T
        c::Vector{T}

        function SolovevEquilibrium{T}(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, a::T, c::Vector{T}) where T <: Number
            new("SolovevEquilibrium", R₀, B₀, ϵ, κ, δ, a, c)
        end
    end

    function SolovevEquilibrium(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, a::T) where T <: Number

        n = 7

        x = [Sym("x" * string(i)) for i in 1:3]
        c = [Sym("c" * string(i)) for i in 1:n]

        ψ = (ψ₀(x,a) + c[1] * ψ₁(x)
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

        SolovevEquilibrium{T}(R₀, B₀, ϵ, κ, δ, a, cnum)
    end


    SolovevEquilibriumITER() = SolovevEquilibrium(6.2, 5.3, 0.32, 1.7, 0.33, -0.155)
    SolovevEquilibriumNSTX() = SolovevEquilibrium(0.85, 0.3, 0.78, 2.0, 0.35, 1.0)


    function Base.show(io::IO, equ::SolovevEquilibrium)
        print(io, "SolovevEquilibrium Equilibrium with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  ϵ  = ", equ.ϵ, "\n")
        print(io, "  κ  = ", equ.κ, "\n")
        print(io, "  δ  = ", equ.δ, "\n")
        print(io, "  a  = ", equ.a)
    end


    function ElectromagneticFields.A₃(x::AbstractArray{T,1}, equ::SolovevEquilibrium) where {T <: Number}
        ψ₀(x, equ.a) + equ.c[1] * ψ₁(x)
                     + equ.c[2] * ψ₂(x)
                     + equ.c[3] * ψ₃(x)
                     + equ.c[4] * ψ₄(x)
                     + equ.c[5] * ψ₅(x)
                     + equ.c[6] * ψ₆(x)
                     + equ.c[7] * ψ₇(x)
    end


    macro solovev_equilibrium(R₀, B₀, ϵ, κ, δ, a)
        generate_equilibrium_code(SolovevEquilibrium(R₀, B₀, ϵ, κ, δ, a); output=false)
    end

    function init(R₀, B₀, ϵ, κ, δ, a; perturbation=ZeroPerturbation())
        equilibrium = SolovevEquilibrium(R₀, B₀, ϵ, κ, δ, a)
        load_equilibrium(equilibrium, perturbation; target_module=Solovev)
        return equilibrium
    end

end
