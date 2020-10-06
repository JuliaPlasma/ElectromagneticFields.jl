@doc raw"""
Axisymmetric Solov'ev equilibra in (R/R₀,Z/R₀,ϕ) coordinates.
Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010,
      and Freidberg, Ideal Magnetohydrodynamics, 2014.

The covariant components of the vector potential are given by
```math
A (x, y, \phi) = \frac{B_0 R_0}{2} \, \left( \frac{y}{x} , \, - \ln x , \, \psi(x,y) \right)^T ,
```

with $x = R/R_0$ and $y = Z/R_0$.
The normalised poloidal flux $\psi$ is given by
```math
\psi (x,y) = \psi_0 + \sum \limits_{i=1}^{7} c_i \psi_i (x,y) ,
```

with
```math
\begin{aligned}
\psi_{0} &= \frac{x^4}{8} + \alpha \left( \frac{1}{2} x^2 \, \ln x - \frac{x^4}{8} \right) , \\
\psi_{1} &= 1 , \\
\psi_{2} &= x^2 , \\
\psi_{3} &= y^2 - x^2 \, \ln x , \\
\psi_{4} &= x^4 - 4 x^2 y^2 , \\
\psi_{5} &= 2 y^4 9 y^2 x^2 + 3 x^4 \, \ln x - 12 x^2 y^2 \, \ln x , \\
\psi_{6} &= x^6 - 12 x^4 y^2 + 8 x^2 y^4 , \\
\psi_{7} &= 8 y^6 - 140 y^4 x^2 + 75 y^2 x^4 - 15 x^6 \, \ln x + 180 x^4 y^2 \, \ln x - 120 x^2 y^4 \, \ln x .
\end{aligned}
```

This formula describes exact solutions of the Grad-Shafranov equation with up-down symmetry.
The constants $c_i$ are determined from boundary constraints on $\psi$, that are derived from
the following analytic model for a smooth, elongated "D" shaped cross section:
```math
\begin{aligned}
x &= 1 + \epsilon \, \cos (\tau + \delta_0 \, \sin \tau) , \\
y &= \epsilon \kappa \, \sin (\tau) ,
\end{aligned}
```
where $0 \leq \tau < 2 \pi$, $\epsilon = a / R_0$ is the inverse aspect ratio, $\kappa$ the elongation,
and $\sin \delta_0 = \delta$ is the triangularity.

Defining three test points, namely
- the high point $(1 - \delta \epsilon, \kappa \epsilon)$,
- the inner equatorial point $(1 - \epsilon, 0)$,
- and the outer equatorial point $(1 + \epsilon, 0)$,

the following geometric constraints can be posed on the solution:
```math
\begin{aligned}
\psi (1 + \epsilon, 0) &= 0 , \\
\psi (1 - \epsilon, 0) &= 0 , \\
\psi (1 - \delta \epsilon, \kappa \epsilon) &= 0 , \\
\psi_{x} (1 - \delta \epsilon, \kappa \epsilon) &= 0 , \\
\psi_{yy} (1 + \epsilon, 0) &= - N_1 \psi_{x} (1 + \epsilon, 0) , \\
\psi_{yy} (1 - \epsilon, 0) &= - N_2 \psi_{x} (1 - \epsilon, 0) , \\
\psi_{xx} (1 - \delta \epsilon, \kappa \epsilon) &= - N_3 \psi_y (1 - \delta \epsilon, \kappa \epsilon) .
\end{aligned}
```

The first three equations define the three test points, the fourth equations enforces the high
point to be a maximum, and the last three equations define the curvature at the test points.

The coefficients $N_j$ can be found from the analytic model cross section as
```math
\begin{aligned}
N_1 &= \left[ \frac{d^2 x}{dy^2} \right]_{\tau = 0}     = - \frac{(1 + \delta_0)^2}{\epsilon \kappa^2} , \\
N_2 &= \left[ \frac{d^2 x}{dy^2} \right]_{\tau = \pi}   = \hphantom{-} \frac{(1 - \delta_0)^2}{\epsilon \kappa^2} , \\
N_3 &= \left[ \frac{d^2 x}{dy^2} \right]_{\tau = \pi/2} = - \frac{\kappa}{\epsilon \, \cos^2 \delta_0} .
\end{aligned}
```

For a given value of the constant $a$ above conditions reduce to a set of seven linear
inhomogeneous algebraic equations for the unknown $c_i$, which can easily be solved.


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


    SolovevEquilibriumITER() = SolovevEquilibrium(6.2, 5.3, 0.32,  1.7, 0.33, -0.155)
    # SolovevEquilibriumTFTR() = SolovevEquilibrium(2.5, 5.6, 0.345, 1.0, 0.0,  )
    # SolovevEquilibriumJET()  = SolovevEquilibrium(3.0, 3.6, 0.333, 1.7, 0.25, )
    SolovevEquilibriumNSTX() = SolovevEquilibrium(0.85, 0.30, 0.78, 2.00, 0.35, 1.0)
    # SolovevEquilibriumMAST() = SolovevEquilibrium(0.85, 0.52, 0.77, 2.45, 0.50,  )
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
