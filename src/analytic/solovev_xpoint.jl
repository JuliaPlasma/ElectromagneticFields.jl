@doc raw"""
Axisymmetric Solov'ev equilibra with X-point in (R/R₀,Z/R₀,phi) coordinates.
Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010,
      and Freidberg, Ideal Magnetohydrodynamics, 2014.

The covariant components of the vector potential are given by
```math
A (x, y, \phi) = \left( \frac{B_0 R_0}{2} \, \frac{y}{x} , \, - \frac{B_0 R_0}{2} \, \ln x , \, \psi(x,y) \right)^T ,
```

with $x = R/R_0$ and $y = Z/R_0$.
The normalised poloidal flux $\psi$ is given by
```math
\psi (x,y) = \psi_0 + \sum \limits_{i=1}^{12} c_i \psi_i (x,y) ,
```

with
```math
\begin{aligned}
\psi_{0}  &= \frac{x^4}{8} + \alpha \left( \frac{1}{2} x^2 \, \ln x - \frac{x^4}{8} \right) , \\
\psi_{1}  &= 1 , \\
\psi_{2}  &= x^2 , \\
\psi_{3}  &= y^2 - x^2 \, \ln x , \\
\psi_{4}  &= x^4 - 4 x^2 y^2 , \\
\psi_{5}  &= 2 y^4 9 y^2 x^2 + 3 x^4 \, \ln x - 12 x^2 y^2 \, \ln x , \\
\psi_{6}  &= x^6 - 12 x^4 y^2 + 8 x^2 y^4 , \\
\psi_{7}  &= 8 y^6 - 140 y^4 x^2 + 75 y^2 x^4 - 15 x^6 \, \ln x + 180 x^4 y^2 \, \ln x - 120 x^2 y^4 \, \ln x , \\
\psi_{8}  &= y , \\
\psi_{9}  &= y x^2 , \\
\psi_{10} &= y^3 - 3 y x^2 \, \ln x , \\
\psi_{11} &= 3 y x^4 - 4 y^3 x^2 , \\
\psi_{12} &= 8 y^5 - 45 y x^4 - 80 y^3 x^2 \, \ln x + 60 y x^4 \, \ln x .
\end{aligned}
```

This formula describes exact solutions of the Grad-Shafranov equation with up-down asymmetry.
The constants $c_i$ are determined from boundary constraints on $\psi$, that are derived from
the following analytic model for a smooth, elongated "D" shaped cross section:
```math
\begin{aligned}
x &= 1 + \epsilon \, \cos (\tau + \arcsin \delta \, \sin \tau) , \\
y &= \epsilon \kappa \, \sin (\tau) ,
\end{aligned}
```
where $0 \leq \tau < 2 \pi$, $\epsilon = a / R_0$ is the inverse aspect ratio, $\kappa$ the elongation,
and $\sin \delta_0 = \delta$ is the triangularity.

Defining four test points, namely
- the high point $(1 - \delta \epsilon, \kappa \epsilon)$,
- the inner equatorial point $(1 - \epsilon, 0)$,
- and the outer equatorial point $(1 + \epsilon, 0)$,
- the position of the X-point $(x_{\mathrm{sep}}, y_{\mathrm{sep}})$,

the following geometric constraints can be posed on the solution:
```math
\begin{aligned}
\psi (1 + \epsilon, 0) &= 0 , \\
\psi (1 - \epsilon, 0) &= 0 , \\
\psi (1 - \delta \epsilon, \kappa \epsilon) &= 0 , \\
\psi (x_{\mathrm{sep}}, y_{\mathrm{sep}}) &= 0 , \\
\psi_{y} (1 + \epsilon, 0) &= 0 , \\
\psi_{y} (1 - \epsilon, 0) &= 0 , \\
\psi_{x} (1 - \delta \epsilon, \kappa \epsilon) &= 0 , \\
\psi_{x} (x_{\mathrm{sep}}, y_{\mathrm{sep}}) &= 0 , \\
\psi_{y} (x_{\mathrm{sep}}, y_{\mathrm{sep}}) &= 0 , \\
\psi_{yy} (1 + \epsilon, 0) &= - N_1 \psi_{x} (1 + \epsilon, 0) , \\
\psi_{yy} (1 - \epsilon, 0) &= - N_2 \psi_{x} (1 - \epsilon, 0) , \\
\psi_{xx} (1 - \delta \epsilon, \kappa \epsilon) &= - N_3 \psi_y (1 - \delta \epsilon, \kappa \epsilon) .
\end{aligned}
```

The first four equations define the four test points, the fifth and sixth equations define the up-down
symmetry, the seventh equations enforces the high point to be a maximum, the eighth and ninth eqaution
set the $x$- and $y$-components of the magnetic field at the X-point to zero, and the last three
equations define the curvature at the first three test points.

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
 * `xsep`: x position of the X point
 * `ysep`: y position of the X point
"""
module SolovevXpoint

    using RecipesBase

    using SymPy: N, Sym, diff, expand, solve, subs

    import ..ElectromagneticFields
    import ..ElectromagneticFields: ZeroPerturbation
    import ..ElectromagneticFields: load_equilibrium, generate_equilibrium_code
    import ..SolovevAbstract: AbstractSolovevEquilibrium, X, Y, Z, R, r, θ, ϕ, r²

    export  SolovevXpointEquilibrium

    include("solovev_psi.jl")

    struct SolovevXpointEquilibrium{T <: Number} <: AbstractSolovevEquilibrium
        name::String
        R₀::T
        B₀::T
        ϵ::T
        κ::T
        δ::T
        α::T
        xsep::T
        ysep::T
        c::Vector{T}

        function SolovevXpointEquilibrium{T}(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, α::T, xsep::T, ysep::T, c::Vector{T}) where T <: Number
            new("Solovev Equilibrium with X-point", R₀, B₀, ϵ, κ, δ, α, xsep, ysep, c)
        end
    end

    function SolovevXpointEquilibrium(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, α::T, xsep::T, ysep::T) where T <: Number

        n = 12

        x = [Sym("x" * string(i)) for i in 1:3]
        c = [Sym("c" * string(i)) for i in 1:n]

        ψ = ( ψ₀(x,α) + c[1]  * ψ₁(x)
                      + c[2]  * ψ₂(x)
                      + c[3]  * ψ₃(x)
                      + c[4]  * ψ₄(x)
                      + c[5]  * ψ₅(x)
                      + c[6]  * ψ₆(x)
                      + c[7]  * ψ₇(x)
                      + c[8]  * ψ₈(x)
                      + c[9]  * ψ₉(x)
                      + c[10] * ψ₁₀(x)
                      + c[11] * ψ₁₁(x)
                      + c[12] * ψ₁₂(x) )

        eqs = [
            expand(subs(subs(ψ, x[1], 1+ϵ), x[2], 0)),
            expand(subs(subs(ψ, x[1], 1-ϵ), x[2], 0)),
            expand(subs(subs(ψ, x[1], 1-δ*ϵ), x[2], κ*ϵ)),
            expand(subs(subs(ψ, x[1], xsep), x[2], ysep)),
            expand(subs(subs(diff(ψ, x[2]), x[1], 1+ϵ), x[2], 0)),
            expand(subs(subs(diff(ψ, x[2]), x[1], 1-ϵ), x[2], 0)),
            expand(subs(subs(diff(ψ, x[1]), x[1], 1-δ*ϵ), x[2], κ*ϵ)),
            expand(subs(subs(diff(ψ, x[1]), x[1], xsep), x[2], ysep)),
            expand(subs(subs(diff(ψ, x[2]), x[1], xsep), x[2], ysep)),
            expand(subs(subs(diff(ψ, x[2], 2), x[1], 1+ϵ), x[2], 0) - (1 + asin(δ))^2 / (ϵ * κ^2) * subs(subs(diff(ψ, x[1]), x[1], 1+ϵ), x[2], 0)),
            expand(subs(subs(diff(ψ, x[2], 2), x[1], 1-ϵ), x[2], 0) + (1 - asin(δ))^2 / (ϵ * κ^2) * subs(subs(diff(ψ, x[1]), x[1], 1-ϵ), x[2], 0)),
            expand(subs(subs(diff(ψ, x[1], 2), x[1], 1-δ*ϵ), x[2], κ*ϵ) - κ / (ϵ * (1 - δ^2)) * subs(subs(diff(ψ, x[2]), x[1], 1-δ*ϵ), x[2], κ*ϵ))
        ]

        csym = solve(eqs, c)
        cnum = [N(csym[c[i]]) for i in 1:n]

        SolovevXpointEquilibrium{T}(R₀, B₀, ϵ, κ, δ, α, xsep, ysep, cnum)
    end


    SolovevXpointEquilibriumITER() = SolovevXpointEquilibrium(6.2, 5.3, 0.32, 1.7, 0.33, -0.155, 0.88, -0.60)
    SolovevXpointEquilibriumNSTX() = SolovevXpointEquilibrium(0.85, 0.3, 0.78, 2.0, 0.35, -0.05, 0.70, -1.71)


    function Base.show(io::IO, equ::SolovevXpointEquilibrium)
        print(io, "Solovev Xpoint Equilibrium with\n")
        print(io, "  R₀ = ", equ.R₀, "\n")
        print(io, "  B₀ = ", equ.B₀, "\n")
        print(io, "  ϵ  = ", equ.ϵ,  "\n")
        print(io, "  κ  = ", equ.κ,  "\n")
        print(io, "  δ  = ", equ.δ,  "\n")
        print(io, "  α  = ", equ.α,  "\n")
        print(io, "  xsep  = ", equ.xsep, "\n")
        print(io, "  ysep  = ", equ.ysep)
    end


    function ElectromagneticFields.A₃(x::AbstractArray{T,1}, equ::SolovevXpointEquilibrium) where {T <: Number}
        ( ψ₀(x, equ.α) + equ.c[1]  * ψ₁(x)
                       + equ.c[2]  * ψ₂(x)
                       + equ.c[3]  * ψ₃(x)
                       + equ.c[4]  * ψ₄(x)
                       + equ.c[5]  * ψ₅(x)
                       + equ.c[6]  * ψ₆(x)
                       + equ.c[7]  * ψ₇(x)
                       + equ.c[8]  * ψ₈(x)
                       + equ.c[9]  * ψ₉(x)
                       + equ.c[10] * ψ₁₀(x)
                       + equ.c[11] * ψ₁₁(x)
                       + equ.c[12] * ψ₁₂(x) )
    end


    macro solovev_xpoint_equilibrium(R₀, B₀, ϵ, κ, δ, α, xsep, ysep)
        generate_equilibrium_code(SolovevXpointEquilibrium(R₀, B₀, ϵ, κ, δ, α, xsep, ysep); output=false)
    end

    function init(R₀, B₀, ϵ, κ, δ, α, xsep, ysep; perturbation=ZeroPerturbation(), target_module=SolovevXpoint)
        equilibrium = SolovevXpointEquilibrium(R₀, B₀, ϵ, κ, δ, α, xsep, ysep)
        load_equilibrium(equilibrium, perturbation; target_module=target_module)
        return equilibrium
    end

    function ITER(; perturbation=ZeroPerturbation(), target_module=SolovevXpoint)
        equilibrium = SolovevXpointEquilibriumITER()
        load_equilibrium(equilibrium, perturbation; target_module=target_module)
        return equilibrium
    end

    function NSTX(; perturbation=ZeroPerturbation(), target_module=SolovevXpoint)
        equilibrium = SolovevXpointEquilibriumNSTX()
        load_equilibrium(equilibrium, perturbation; target_module=target_module)
        return equilibrium
    end


    @recipe function f(equ::SolovevXpointEquilibrium;
                       nx = 100, ny = 120, levels = 50, size = (300,400),
                       xlims = ( 0.50,  1.50),
                       ylims = (-0.75, +0.75))

        xgrid = LinRange(xlims[1], xlims[2], nx)
        zgrid = LinRange(ylims[1], ylims[2], ny)
        pot   = [A₃(0, xgrid[i], zgrid[j], 0.0) / xgrid[i] for i in eachindex(xgrid), j in eachindex(zgrid)]

        seriestype := :contour
        aspect_ratio := :equal
        size   := size
        xlims  := xlims
        ylims  := ylims
        levels := levels
        legend := :none

        (xgrid, zgrid, pot')
    end
    
end
