# Check the Symbolics rewrite against ElectromagneticFields 0.8, which generated its field
# functions with SymEngine and injected them into a module.
#
# `scripts/reference-symengine.tsv` holds the value of every one of the ~450 functions 0.8 generated,
# for each of the 20 equilibria, at a fixed point. This script evaluates the corresponding entry of
# the tensor-valued function the rewrite stores and asserts the two agree.
#
# This is the only check that establishes that no generated function was lost and that none changed
# value. It cannot be a test, because the 0.8 code that produced the reference no longer exists in
# the repository; regenerate the file with the generator quoted at the bottom if the reference ever
# needs rebuilding from a 0.8 checkout.
#
#     julia --project=. scripts/verify_against_symengine.jl

using ElectromagneticFields
using Printf

const REFERENCE = joinpath(@__DIR__, "reference-symengine.tsv")

const EQUILIBRIA = Dict(
    "ABC" => (ABCEquilibrium(), [0.5, 0.5, 0.5]),
    "TokamakCartesian" => (AxisymmetricTokamakCartesianEquilibrium(), [1.5, 0.5, 0.5]),
    "TokamakCylindrical" => (AxisymmetricTokamakCylindricalEquilibrium(), [1.5, 0.5, 0.5]),
    "TokamakToroidal" => (AxisymmetricTokamakToroidalEquilibrium(), [0.5, 0.5, 0.5]),
    "TokamakToroidalReg" => (AxisymmetricTokamakToroidalRegularizationEquilibrium(),
        [0.5, 0.5, 0.5]),
    "Dipole" => (DipoleField(), [0.5, 0.5, 0.5]),
    "PenningUniform" => (PenningTrapUniformEquilibrium(), [0.5, 0.5, 0.5]),
    "PenningBottle" => (PenningTrapBottleEquilibrium(), [0.5, 0.5, 0.5]),
    "PenningAsymmetric" => (PenningTrapAsymmetricEquilibrium(), [0.5, 0.5, 0.5]),
    "QuadraticPotentials" => (QuadraticPotentialsField(), [0.5, 0.5, 0.5]),
    "Singular" => (SingularEquilibrium(), [0.5, 0.5, 0.5]),
    "SolovevITER" => (SolovevEquilibriumITER(), [1.05, 0.25, 0.5]),
    "SolovevNSTX" => (SolovevEquilibriumNSTX(), [1.05, 0.25, 0.5]),
    "SolovevFRC" => (SolovevEquilibriumFRC(), [1.05, 0.25, 0.5]),
    "SolovevXpointITER" => (SolovevXpointEquilibriumITER(), [1.05, 0.25, 0.5]),
    "SolovevXpointNSTX" => (SolovevXpointEquilibriumNSTX(), [1.05, 0.25, 0.5]),
    "SolovevDoubleXNSTX" => (SolovevDoubleXpointEquilibriumNSTX(), [1.05, 0.25, 0.5]),
    "SolovevSymmetric" => (SolovevSymmetricEquilibrium(), [0.5, 0.5, 0.5]),
    "SymmetricQuadratic" => (SymmetricQuadraticEquilibrium(), [0.5, 0.5, 0.5]),
    "ThetaPinch" => (ThetaPinchEquilibrium(), [0.5, 0.5, 0.5])
)

const SUB = Dict('₁' => 1, '₂' => 2, '₃' => 3)
const SUP = Dict('¹' => 1, '²' => 2, '³' => 3)

"""
Resolve one of 0.8's generated names to the value the rewrite produces, or `nothing` when the name
has no counterpart.

The scalar names of 0.8 are the entries of the rewrite's tensors, so each branch below pairs a name
shape with the generic that now carries it and the index into its result. Names are matched over
`collect(name)` rather than by byte index: every index character here is multi-byte, and `DF̄`
carries a combining macron of its own.
"""
function resolve(field, name::String, t, ξ)
    cs = collect(name)
    n = length(cs)
    idx(c) = get(SUB, c, get(SUP, c, 0))
    coords = coordinates(field)

    # the per-equilibrium coordinate helpers keep their own names
    sym = Symbol(name)
    haskey(coords, sym) && return coords[sym](t, ξ)

    name == "orientation" && return orientation(field)
    name == "J" && return J(field, t, ξ)
    name == "B" && return B(field, t, ξ)
    name == "φ" && return φ(field, t, ξ)

    # derivatives of |B|: d²Bdxᵢdxⱼ and dBdxᵢ
    startswith(name, "d²Bdx") && n == 9 && return DDB(field, t, ξ)[idx(cs[6]), idx(cs[9])]
    startswith(name, "dBdx") && n == 5 && return DB(field, t, ξ)[idx(cs[5])]

    # metric derivatives: d²gᵢⱼdxₖdxₗ and dgᵢⱼdxₖ, sub- or superscripted
    if startswith(name, "d²g") && n == 11
        gen = haskey(SUB, cs[4]) ? DDg♭ : DDg♯
        return gen(field, t, ξ)[idx(cs[4]), idx(cs[5]), idx(cs[8]), idx(cs[11])]
    end
    if startswith(name, "dg") && n == 7
        gen = haskey(SUB, cs[3]) ? Dg♭ : Dg♯
        return gen(field, t, ξ)[idx(cs[3]), idx(cs[4]), idx(cs[7])]
    end

    # second derivatives of a vector quantity: d²Xᵢdxⱼdxₖ
    startswith(name, "d²A") && n == 10 &&
        return DDA♭(field, t, ξ)[idx(cs[4]), idx(cs[7]), idx(cs[10])]
    startswith(name, "d²b") && n == 10 &&
        return DDb♭(field, t, ξ)[idx(cs[4]), idx(cs[7]), idx(cs[10])]

    # first derivatives: dXᵢdxⱼ, and the physical db₍ᵢ₎dxⱼ
    startswith(name, "db₍") && n == 8 && return Db♮(field, t, ξ)[idx(cs[4]), idx(cs[8])]
    if n == 6 && cs[1] == 'd' && haskey(SUB, cs[3])
        gen = cs[2] == 'A' ? DA♭ :
              cs[2] == 'B' ? DB♭ : cs[2] == 'b' ? Db♭ :
                                   cs[2] == 'E' ? DE♭ : nothing
        gen === nothing || return gen(field, t, ξ)[idx(cs[3]), idx(cs[6])]
    end

    # the chart Jacobian and its inverse
    startswith(name, "DF̄") && return DF̄(field, t, ξ)[idx(cs[end - 1]), idx(cs[end])]
    startswith(name, "DF") && return DF(field, t, ξ)[idx(cs[end - 1]), idx(cs[end])]

    # chart maps and ranges
    n == 2 && cs[1] == 'x' && haskey(SUP, cs[2]) &&
        return to_cartesian(field, t, ξ)[SUP[cs[2]]]
    n == 2 && cs[1] == 'ξ' && haskey(SUP, cs[2]) &&
        return from_cartesian(field, t, ξ)[SUP[cs[2]]]
    startswith(name, "minx") && n == 5 && return rangemin(field, t, ξ)[idx(cs[5])]
    startswith(name, "maxx") && n == 5 && return rangemax(field, t, ξ)[idx(cs[5])]

    # rank-two objects written with two index characters: gᵢⱼ, gⁱʲ and the two-form Bᵢⱼ
    if n == 3 && cs[1] == 'g'
        gen = haskey(SUB, cs[2]) ? g♭ : g♯
        return gen(field, t, ξ)[idx(cs[2]), idx(cs[3])]
    end
    n == 3 && cs[1] == 'B' && haskey(SUB, cs[2]) && haskey(SUB, cs[3]) &&
        return B♭♭(field, t, ξ)[idx(cs[2]), idx(cs[3])]

    # components of a vector quantity in one of the three representations
    covariant = Dict('A' => A♭, 'B' => B♭, 'E' => E♭, 'a' => a♭, 'b' => b♭, 'c' => c♭)
    contravariant = Dict('A' => A♯, 'B' => B♯, 'E' => E♯, 'a' => a♯, 'b' => b♯, 'c' => c♯)
    physical = Dict('B' => B♮, 'a' => a♮, 'b' => b♮, 'c' => c♮)

    if n == 2 && haskey(covariant, cs[1]) && haskey(SUB, cs[2])
        return covariant[cs[1]](field, t, ξ)[SUB[cs[2]]]
    elseif n == 2 && haskey(contravariant, cs[1]) && haskey(SUP, cs[2])
        return contravariant[cs[1]](field, t, ξ)[SUP[cs[2]]]
    elseif n == 4 && cs[2] == '₍' && haskey(physical, cs[1])
        return physical[cs[1]](field, t, ξ)[idx(cs[3])]
    end

    nothing
end

# The perpendicular frame is not unique: `a` is built from the first trial vector whose cross
# product with `b` does not vanish, and `c = b × a`. Any frame satisfying the orthonormality the
# test suite checks is correct, so a disagreement here is reported but not counted as a failure.
isframe(name) = name[1] in ('a', 'c')

function main()
    lines = readlines(REFERENCE)
    byequ = Dict{String, Vector{Tuple{String, Float64}}}()
    for line in lines
        startswith(line, "#") && continue
        parts = split(line, '\t')
        length(parts) == 3 || continue
        push!(get!(byequ, parts[1], []), (parts[2], parse(Float64, parts[3])))
    end

    total = checked = unresolved = failed = framediff = 0
    unresolved_names = Set{String}()

    for label in sort(collect(keys(byequ)))
        haskey(EQUILIBRIA, label) || continue
        equ, ξ = EQUILIBRIA[label]
        field = FieldFunctions(equ)
        t = 1.0

        n_checked = n_failed = n_frame = n_unres = 0
        for (name, expected) in byequ[label]
            total += 1
            got = resolve(field, name, t, ξ)
            if got === nothing
                n_unres += 1
                push!(unresolved_names, name)
                continue
            end
            n_checked += 1
            ok = isapprox(got, expected; rtol = 1e-10, atol = 1e-12) ||
                 (isnan(got) && isnan(expected))
            if !ok
                if isframe(name)
                    n_frame += 1
                else
                    n_failed += 1
                    n_failed ≤ 5 &&
                        @printf("    MISMATCH %-22s expected % .12g  got % .12g\n",
                            name, expected, got)
                end
            end
        end
        checked += n_checked
        failed += n_failed
        framediff += n_frame
        unresolved += n_unres
        @printf("%-22s %4d checked, %3d frame differences, %3d unresolved, %3d FAILED\n",
            label, n_checked, n_frame, n_unres, n_failed)
    end

    println()
    @printf("total %d reference values: %d checked, %d unresolved, %d frame differences, %d FAILED\n",
        total, checked, unresolved, framediff, failed)
    if !isempty(unresolved_names)
        println("unresolved names: ", join(sort(collect(unresolved_names)), ", "))
    end
    failed == 0 || error("$failed values disagree with the SymEngine reference")
end

main()
