
using SymPy: N, Sym, diff, expand, solve, subs

"""
Axisymmetric Solov'ev equilibra with X-point in (R,Z,phi) coordinates.
Based on Cerfon & Freidberg, Physics of Plasmas 17, 032502, 2010.

Parameters:
    R₀: position of magnetic axis
    B₀: B-field at magnetic axis
    ϵ:  inverse aspect ratio
    κ:  elongation
    δ:  triangularity
    a:  free constant, determined to match a given beta value
    xₛₑₚ: x position of the X point
    yₛₑₚ: y position of the X point
"""
struct SolovevXpoint{T <: Number} <: AnalyticEquilibrium
    const name::String = "SolovevXpointEquilibrium"
    R₀::T
    B₀::T
    ϵ::T
    κ::T
    δ::T
    a::T
    xₛₑₚ::T
    yₛₑₚ::T
    c::Vector{T}

    function SolovevXpoint{T}(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, a::T, xₛₑₚ::T, yₛₑₚ::T, c::Vector{T}) where T <: Number
        new(R₀, B₀, ϵ, κ, δ, a, xₛₑₚ, yₛₑₚ, c)
    end
end

function SolovevXpoint(R₀::T, B₀::T, ϵ::T, κ::T, δ::T, a::T, xₛₑₚ::T, yₛₑₚ::T) where T <: Number

    n = 12

    x = [Sym("x" * string(i)) for i in 1:3]
    c = [Sym("c" * string(i)) for i in 1:n]

    ψ = (ψ₀(x,a) + c[1]  * ψ₁(x)
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
        expand(subs(subs(ψ, x[1], xₛₑₚ), x[2], yₛₑₚ)),
        expand(subs(subs(diff(ψ, x[2]), x[1], 1+ϵ), x[2], 0)),
        expand(subs(subs(diff(ψ, x[2]), x[1], 1-ϵ), x[2], 0)),
        expand(subs(subs(diff(ψ, x[1]), x[1], 1-δ*ϵ), x[2], κ*ϵ)),
        expand(subs(subs(diff(ψ, x[1]), x[1], xₛₑₚ), x[2], yₛₑₚ)),
        expand(subs(subs(diff(ψ, x[2]), x[1], xₛₑₚ), x[2], yₛₑₚ)),
        expand(subs(subs(diff(ψ, x[2], 2), x[1], 1+ϵ), x[2], 0) - (1 + asin(δ))^2 / (ϵ * κ^2) * subs(subs(diff(ψ, x[1]), x[1], 1+ϵ), x[2], 0)),
        expand(subs(subs(diff(ψ, x[2], 2), x[1], 1-ϵ), x[2], 0) + (1 - asin(δ))^2 / (ϵ * κ^2) * subs(subs(diff(ψ, x[1]), x[1], 1-ϵ), x[2], 0)),
        expand(subs(subs(diff(ψ, x[1], 2), x[1], 1-δ*ϵ), x[2], κ*ϵ) - κ / (ϵ * (1 - δ^2)) * subs(subs(diff(ψ, x[2]), x[1], 1-δ*ϵ), x[2], κ*ϵ))
    ]

    csym = solve(eqs, c)
    cnum = [N(csym[c[i]]) for i in 1:n]

    SolovevXpoint{T}(R₀, B₀, ϵ, κ, δ, a, xₛₑₚ, yₛₑₚ, cnum)
end


function Base.show(io::IO, equ::SolovevXpoint)
    print(io, "Solovev Xpoint Equilibrium with\n")
    print(io, "  R₀ = ", equ.R₀, "\n")
    print(io, "  B₀ = ", equ.B₀, "\n")
    print(io, "  ϵ  = ", equ.ϵ, "\n")
    print(io, "  κ  = ", equ.κ, "\n")
    print(io, "  δ  = ", equ.δ, "\n")
    print(io, "  a  = ", equ.a)
end


function analyticA₁(x::Vector, equ::SolovevXpoint)
    + equ.B₀ * equ.R₀ * x[2] / x[1] / 2
end

function analyticA₂(x::Vector, equ::SolovevXpoint)
    - equ.B₀ * equ.R₀ * log(x[1]) / 2
end

function analyticA₃(x::Vector, equ::SolovevXpoint)
    (ψ₀(x, equ.a) + equ.c[1]  * ψ₁(x)
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

function analyticMetric(x::Vector, equ::SolovevXpoint)
    R = x[1] * equ.R₀
    [1  0  0;
     0  1  0;
     0  0  R^2]
end


macro solovev_xpoint_equilibrium(R₀, B₀, ϵ, κ, δ, a, xₛₑₚ, yₛₑₚ)
    generate_equilibrium_code(SolovevXpoint(R₀, B₀, ϵ, κ, δ, a, xₛₑₚ, yₛₑₚ); output=false)
end
