# Check the Penning trap docstrings against the code they document.
#
# Each of the three Penning trap docstrings states the vector potential `A`, the magnetic field
# `B = ∇ × A`, the electrostatic potential `φ` and the electric field `E = -∇φ` in closed form.
# This script rebuilds `A` and `φ` from the package's own methods with symbolic parameters, takes
# the curl and the gradient, and compares all four against the formulas the docstrings give.
#
# The comparison is numeric, at random points: the docstring formulas are transcribed with `/`
# as the docstrings write them, so a symbolic difference can carry a `Float64` coefficient where
# the code carries a `Rational` and not reduce to a structural zero.
#
# Run it from the repository root:
#
#     julia --startup-file=no --project=. scripts/verify_penning_trap_docstrings.jl

using ElectromagneticFields
using ElectromagneticFields: A₁, A₂, A₃, φ
using Random
using Symbolics

Symbolics.@variables x y z B₀ Bₚ E₀

const ξ = [x, y, z]
const D = [Symbolics.Differential(v) for v in ξ]
const VARIABLES = [x, y, z, B₀, Bₚ, E₀]
const NPOINTS = 8
const RTOL = 1e-10

curl(A) = [D[2](A[3]) - D[3](A[2]),
    D[3](A[1]) - D[1](A[3]),
    D[1](A[2]) - D[2](A[1])]

grad(f) = [D[i](f) for i in 1:3]

# `expand_derivatives` resolves the `Differential` terms the curl and the gradient introduce.
# `simplify` is deliberately not used, for the reason given at the `iszero` guard in
# `src/analytic/analytic_field.jl`.
resolved(u) = Symbolics.expand_derivatives.(u)

"""
Whether `u` and `v` agree, judged at `NPOINTS` random points, relative to the larger of the two.
"""
function agrees(rng, u, v)
    du, dv = resolved(u), resolved(v)

    for _ in 1:NPOINTS
        point = Dict(VARIABLES .=> randn(rng, length(VARIABLES)))
        a = Symbolics.value.(Symbolics.substitute.(du, Ref(point)))
        b = Symbolics.value.(Symbolics.substitute.(dv, Ref(point)))
        scale = max(maximum(abs, a), maximum(abs, b), one(eltype(a)))
        all(abs.(a .- b) .<= RTOL * scale) || return false
    end

    true
end

function check(rng, name, equ, A_doc, B_doc, φ_doc, E_doc)
    A_code = [A₁(ξ, equ), A₂(ξ, equ), A₃(ξ, equ)]
    φ_code = φ(ξ, equ)

    results = ["A" => agrees(rng, A_code, A_doc),
        "B" => agrees(rng, curl(A_code), B_doc),
        "φ" => agrees(rng, [φ_code], [φ_doc]),
        "E" => agrees(rng, -grad(φ_code), E_doc)]

    for (quantity, ok) in results
        println(rpad(name, 24), rpad(quantity, 4), ok ? "agrees" : "DISAGREES")
    end

    all(last, results)
end

rng = Random.MersenneTwister(0)

# φ and E are stated identically in all three docstrings.
φ_doc = -E₀ * (x^2 / 2 + y^2 / 2 - z^2)
E_doc = E₀ * [x, y, -2z]

ok = check(rng, "PenningTrapUniform", PenningTrapUniformEquilibrium{Num}(B₀, E₀),
    B₀ * [0, x, 0],
    B₀ * [0, 0, 1],
    φ_doc, E_doc)

ok &= check(rng, "PenningTrapBottle", PenningTrapBottleEquilibrium{Num}(B₀, Bₚ, E₀),
    B₀ / 2 * [-y, x, 0] - Bₚ * [y * z^2 - y^3 / 6, x^3 / 6, x * y * z],
    B₀ * [0, 0, 1] - Bₚ * [x * z, y * z, (x^2 + y^2) / 2 - z^2],
    φ_doc, E_doc)

ok &= check(
    rng, "PenningTrapAsymmetric", PenningTrapAsymmetricEquilibrium{Num}(B₀, Bₚ, E₀),
    B₀ / 2 * [-y, x - z / 6, y / 6] + Bₚ / 2 * [z^2 - y^2, z^2 - x^2, y^2 - x^2],
    B₀ * [1 / 6, 0, 1] + Bₚ * [y - z, x + z, y - x],
    φ_doc, E_doc)

println()
println(ok ? "all docstring formulas agree with the code" : "SOME FORMULAS DISAGREE")
exit(ok ? 0 : 1)
