"""
Axisymmetric Solov'ev equilibra in (R/R₀,Z/R₀,ϕ) coordinates.
Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010,
      and Freidberg, Ideal Magnetohydrodynamics, 2014.
"""
module Solovev

    using RecipesBase

    using SymEngine: N, symbols, diff, expand, subs

    import ..ElectromagneticFields
    import ..ElectromagneticFields: code
    import ..SolovevAbstract: AbstractSolovevEquilibrium, X, Y, Z, R, r, θ, ϕ, r²

    export  SolovevEquilibrium, SolovevXpointEquilibrium


    function ψ₀(x::AbstractVector{T}, a) where {T <: Number}
        x[1]^4 / 8 + a * (x[1]^2 * log(x[1]) / 2 - x[1]^4 / 8 )
    end

    function ψ₁(x::AbstractVector{T}) where {T <: Number}
        one(T)
    end

    function ψ₂(x::AbstractVector{T}) where {T <: Number}
        x[1]^2
    end

    function ψ₃(x::AbstractVector{T}) where {T <: Number}
        x[2]^2 - x[1]^2 * log(x[1])
    end

    function ψ₄(x::AbstractVector{T}) where {T <: Number}
        x[1]^4 - 4 * x[1]^2 * x[2]^2
    end

    function ψ₅(x::AbstractVector{T}) where {T <: Number}
        2 * x[2]^4 - 9 * x[2]^2 * x[1]^2 + 3 * x[1]^4 * log(x[1]) - 12 * x[1]^2 * x[2]^2 * log(x[1])
    end

    function ψ₆(x::AbstractVector{T}) where {T <: Number}
        x[1]^6 - 12 * x[1]^4 * x[2]^2 + 8 * x[1]^2 * x[2]^4
    end

    function ψ₇(x::AbstractVector{T}) where {T <: Number}
        8 * x[2]^6 - 140 * x[2]^4 * x[1]^2 + 75 * x[2]^2 * x[1]^4 - 15 * x[1]^6 * log(x[1]) +
                    180 * x[1]^4 * x[2]^2 * log(x[1]) - 120 * x[1]^2 * x[2]^4 * log(x[1])
    end

    function ψ₈(x::AbstractVector{T}) where {T <: Number}
        x[2]
    end

    function ψ₉(x::AbstractVector{T}) where {T <: Number}
        x[2] * x[1]^2
    end

    function ψ₁₀(x::AbstractVector{T}) where {T <: Number}
        x[2]^3 - 3 * x[2] * x[1]^2 * log(x[1])
    end

    function ψ₁₁(x::AbstractVector{T}) where {T <: Number}
        3 * x[2] * x[1]^4 - 4 * x[2]^3 * x[1]^2
    end

    function ψ₁₂(x::AbstractVector{T}) where {T <: Number}
        8 * x[2]^5 - 45 * x[2] * x[1]^4 - 80 * x[2]^3 * x[1]^2 * log(x[1]) + 60 * x[2] * x[1]^4 * log(x[1])
    end
 

    @doc raw"""
    Axisymmetric Solov'ev equilibra in (R/R₀,Z/R₀,ϕ) coordinates.
    Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010,
        and Freidberg, Ideal Magnetohydrodynamics, 2014.

    The covariant components of the vector potential are given by
    ```math
    A (x, y, \phi) = \left( \frac{B_0 R_0}{2} \, \frac{y}{x} , \, - \frac{B_0 R_0}{2} \, \ln x , \, \psi(x,y) \right)^T ,
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
        A = zeros(n,n+1)

        x₁, x₂, x₃ = symbols("x₁, x₂, x₃")
        x = [x₁, x₂, x₃]
        ψ = [ψ₁(x), ψ₂(x), ψ₃(x), ψ₄(x), ψ₅(x), ψ₆(x), ψ₇(x), -ψ₀(x,α)]

        for i in axes(A,2)
            A[1,i] = N(subs(ψ[i], x[1]=>1+ϵ, x[2]=>0))
            A[2,i] = N(subs(ψ[i], x[1]=>1-ϵ, x[2]=>0))
            A[3,i] = N(subs(ψ[i], x[1]=>1-δ*ϵ, x[2]=>κ*ϵ))
            A[4,i] = N(subs(diff(ψ[i], x[1]), x[1]=>1-δ*ϵ, x[2]=>κ*ϵ))
            A[5,i] = N(subs(diff(ψ[i], x[2], 2), x[1]=>1+ϵ, x[2]=>0) - (1 + asin(δ))^2 / (ϵ * κ^2) * subs(diff(ψ[i], x[1]), x[1]=>1+ϵ, x[2]=>0))
            A[6,i] = N(subs(diff(ψ[i], x[2], 2), x[1]=>1-ϵ, x[2]=>0) + (1 - asin(δ))^2 / (ϵ * κ^2) * subs(diff(ψ[i], x[1]), x[1]=>1-ϵ, x[2]=>0))
            A[7,i] = N(subs(diff(ψ[i], x[1], 2), x[1]=>1-δ*ϵ, x[2]=>κ*ϵ) - κ / (ϵ * (1 - δ^2)) * subs(diff(ψ[i], x[2]), x[1]=>1-δ*ϵ, x[2]=>κ*ϵ))
        end

        c = A[1:n,1:n] \ A[1:n,n+1]

        SolovevEquilibrium{T}(R₀, B₀, ϵ, κ, δ, α, c)
    end

    function init(R₀, B₀, ϵ, κ, δ, α)
        SolovevEquilibrium(R₀, B₀, ϵ, κ, δ, α)
    end

    macro code(R₀, B₀, ϵ, κ, δ, α)
        code(init(R₀, B₀, ϵ, κ, δ, α); escape=true)
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
        A = zeros(n,n+1)

        x₁, x₂, x₃ = symbols("x₁, x₂, x₃")
        x = [x₁, x₂, x₃]
        ψ = [ψ₁(x), ψ₂(x), ψ₃(x), ψ₄(x), ψ₅(x), ψ₆(x), ψ₇(x), ψ₈(x), ψ₉(x), ψ₁₀(x), ψ₁₁(x), ψ₁₂(x), -ψ₀(x,α)]

        for i in axes(A,2)
            A[ 1,i] = N(subs(ψ[i], x[1]=>1+ϵ, x[2]=>0))
            A[ 2,i] = N(subs(ψ[i], x[1]=>1-ϵ, x[2]=>0))
            A[ 3,i] = N(subs(ψ[i], x[1]=>1-δ*ϵ, x[2]=>κ*ϵ))
            A[ 4,i] = N(subs(ψ[i], x[1]=>xsep, x[2]=>ysep))
            A[ 5,i] = N(subs(diff(ψ[i], x[2]), x[1]=>1+ϵ, x[2]=>0))
            A[ 6,i] = N(subs(diff(ψ[i], x[2]), x[1]=>1-ϵ, x[2]=>0))
            A[ 7,i] = N(subs(diff(ψ[i], x[1]), x[1]=>1-δ*ϵ, x[2]=>κ*ϵ))
            A[ 8,i] = N(subs(diff(ψ[i], x[1]), x[1]=>xsep, x[2]=>ysep))
            A[ 9,i] = N(subs(diff(ψ[i], x[2]), x[1]=>xsep, x[2]=>ysep))
            A[10,i] = N(subs(diff(ψ[i], x[2], 2), x[1]=>1+ϵ, x[2]=>0) - (1 + asin(δ))^2 / (ϵ * κ^2) * subs(diff(ψ[i], x[1]), x[1]=>1+ϵ, x[2]=>0))
            A[11,i] = N(subs(diff(ψ[i], x[2], 2), x[1]=>1-ϵ, x[2]=>0) + (1 - asin(δ))^2 / (ϵ * κ^2) * subs(diff(ψ[i], x[1]), x[1]=>1-ϵ, x[2]=>0))
            A[12,i] = N(subs(diff(ψ[i], x[1], 2), x[1]=>1-δ*ϵ, x[2]=>κ*ϵ) - κ / (ϵ * (1 - δ^2)) * subs(diff(ψ[i], x[2]), x[1]=>1-δ*ϵ, x[2]=>κ*ϵ))
        end

        c = A[1:n,1:n] \ A[1:n,n+1]

        SolovevXpointEquilibrium{T}(R₀, B₀, ϵ, κ, δ, α, xsep, ysep, c)
    end


    function SolovevDoubleXpointEquilibrium(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, α::T, xsep::T, ysep::T) where T <: Number

        n = 7
        A = zeros(n,n+1)

        x₁, x₂, x₃ = symbols("x₁, x₂, x₃")
        x = [x₁, x₂, x₃]
        ψ = [ψ₁(x), ψ₂(x), ψ₃(x), ψ₄(x), ψ₅(x), ψ₆(x), ψ₇(x), -ψ₀(x,α)]

        for i in axes(A,2)
            A[1,i] = N(subs(ψ[i], x[1]=>1+ϵ, x[2]=>0))
            A[2,i] = N(subs(ψ[i], x[1]=>1-ϵ, x[2]=>0))
            A[3,i] = N(subs(ψ[i], x[1]=>xsep, x[2]=>ysep))
            A[4,i] = N(subs(diff(ψ[i], x[1]), x[1]=>xsep, x[2]=>ysep))
            A[5,i] = N(subs(diff(ψ[i], x[2]), x[1]=>xsep, x[2]=>ysep))
            A[6,i] = N(subs(diff(ψ[i], x[2], 2), x[1]=>1+ϵ, x[2]=>0) - (1 + asin(δ))^2 / (ϵ * κ^2) * subs(diff(ψ[i], x[1]), x[1]=>1+ϵ, x[2]=>0))
            A[7,i] = N(subs(diff(ψ[i], x[2], 2), x[1]=>1-ϵ, x[2]=>0) + (1 - asin(δ))^2 / (ϵ * κ^2) * subs(diff(ψ[i], x[1]), x[1]=>1-ϵ, x[2]=>0))
        end

        c = vcat(A[1:n,1:n] \ A[1:n,n+1], zeros(5))

        SolovevXpointEquilibrium{T}(R₀, B₀, ϵ, κ, δ, α, xsep, ysep, c)
    end

    function init(R₀, B₀, ϵ, κ, δ, α, xsep, ysep, doublex=false)
        if doublex
            return SolovevDoubleXpointEquilibrium(R₀, B₀, ϵ, κ, δ, α, xsep, ysep)
        else
            return SolovevXpointEquilibrium(R₀, B₀, ϵ, κ, δ, α, xsep, ysep)
        end
    end

    macro code_xpoint(R₀, B₀, ϵ, κ, δ, α, xsep, ysep, doublex=false)
        code(SolovevXpointEquilibrium(R₀, B₀, ϵ, κ, δ, α, xsep, ysep, doublex); escape=true)
    end


    SolovevXpointEquilibriumITER() = SolovevXpointEquilibrium(6.2, 5.3, 0.32, 1.7, 0.33, -0.155, 0.88, -0.60)
    SolovevXpointEquilibriumNSTX() = SolovevXpointEquilibrium(0.85, 0.3, 0.78, 2.0, 0.35, -0.05, 0.70, -1.71)
    SolovevDoubleXpointEquilibriumNSTX() = SolovevDoubleXpointEquilibrium(0.85, 0.3, 0.78, 2.0, 0.35, 0.0, 0.70, -1.71)


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


    function ITER(; xpoint=false)
        if xpoint
            return SolovevXpointEquilibriumITER()
        else
            return SolovevEquilibriumITER()
        end
    end

    function NSTX(; xpoint=false)
        if xpoint
            return SolovevXpointEquilibriumNSTX()
        else
            return SolovevEquilibriumNSTX()
        end
    end

    function NSTXdoubleX()
        SolovevDoubleXpointEquilibriumNSTX()
    end

    function FRC()
        SolovevEquilibriumFRC()
    end


    macro code_iter(xpoint=false)
        code(ITER(xpoint=xpoint); escape=true)
    end

    macro code_iter_xpoint()
        code(ITER(xpoint=true); escape=true)
    end

    macro code_nstx(xpoint=false)
        code(NSTX(xpoint=xpoint); escape=true)
    end

    macro code_nstx_xpoint(doublex=false)
        if doublex
            equilibrium = NSTXdoubleX()
        else
            equilibrium = NSTX(xpoint=true)
        end
        code(equilibrium; escape=true)
    end

    macro code_nstx_double_xpoint()
        code(NSTXdoubleX(); escape=true)
    end

    macro code_frc()
        code(FRC(); escape=true)
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


module SolovevFRC
    import ..Solovev: @code_frc
    export @code
    var"@code" = var"@code_frc"
end

module SolovevITER
    import ..Solovev: @code_iter
    export @code
    var"@code" = var"@code_iter"
end

module SolovevITERwXpoint
    import ..Solovev: @code_iter_xpoint
    export @code
    var"@code" = var"@code_iter_xpoint"
end

module SolovevNSTX
    import ..Solovev: @code_nstx
    export @code
    var"@code" = var"@code_nstx"
end

module SolovevNSTXwXpoint
    import ..Solovev: @code_nstx_xpoint
    export @code
    var"@code" = var"@code_nstx_xpoint"
end

module SolovevNSTXwDoubleXpoint
    import ..Solovev: @code_nstx_double_xpoint
    export @code
    var"@code" = var"@code_nstx_double_xpoint"
end
