# Diagnostic for ChargedParticleDynamics K8, not for merging: after a fresh precompile of this
# package, compare every entry of `DDb♭` with a central difference of `Db♭`, for the fields that
# ChargedParticleDynamics' failing tests use, and save the generated code of `Db♭` and `DDb♭`.
# Usage: julia --project=. probe.jl <output directory>; exit code 2 when an entry is wrong.
using ElectromagneticFields
using ElectromagneticFields: RuntimeGeneratedFunctions

const OUT = ARGS[1]
mkpath(OUT)

function findrgf(x, depth = 0)
    depth > 6 && return nothing
    x isa RuntimeGeneratedFunctions.RuntimeGeneratedFunction && return x
    for n in fieldnames(typeof(x))
        isdefined(x, n) || continue
        r = findrgf(getfield(x, n), depth + 1)
        r === nothing || return r
    end
    nothing
end

const CASES = (
    ("SolovevXpointITER", SolovevXpointEquilibriumITER(),
        ([2.5, 0.0, 0.0], [2.6, 0.1, 0.2])),
    ("TokamakSmallCartesian", AxisymmetricTokamakCartesianEquilibrium(),
        ([1.05, 0.0, 0.0], [1.1, 0.05, 0.02])),
    ("TokamakMediumCartesian", AxisymmetricTokamakCartesianEquilibrium(2.0, 5.0, 2.0),
        ([2.1, 0.0, 0.0], [2.2, 0.1, 0.1])),
    ("TokamakSmallToroidal", AxisymmetricTokamakToroidalEquilibrium(),
        ([1.05, 0.0, 0.0], [1.1, 0.05, 0.02])))

total = 0
for (label, equ, points) in CASES
    field = FieldFunctions(equ)
    for name in (:Db♭, :DDb♭)
        r = findrgf(getfield(field.functions, name))
        body = r === nothing ? "no RuntimeGeneratedFunction found" :
               string(RuntimeGeneratedFunctions.get_expression(r))
        h = string(hash(body); base = 16)
        write(joinpath(OUT, "$(label)_$(name)_$(h).jl"), body)
        println("[k8] ", label, " ", name, " code hash ", h)
    end
    bad = 0
    for x in points
        q = collect(from_cartesian(field, 0.0, x))
        D = DDb♭(field, 0.0, q)
        for i in 1:3, j in 1:3, k in 1:3
            δ = 1e-6 * max(1.0, abs(q[k]))
            qp = copy(q)
            qp[k] += δ
            qm = copy(q)
            qm[k] -= δ
            fd = (Db♭(field, 0.0, qp)[i, j] - Db♭(field, 0.0, qm)[i, j]) / 2δ
            if !isapprox(D[i, j, k], fd; rtol = 1e-5, atol = 1e-8)
                bad += 1
                println("[k8] BAD ", label, " q=", q, " ",
                    (i, j, k), " DDb♭=", D[i, j, k], " fd=", fd)
            end
        end
    end
    println("[k8] ", label, " wrong entries: ", bad)
    global total += bad
end
println("[k8] total wrong entries: ", total)
exit(total == 0 ? 0 : 2)
