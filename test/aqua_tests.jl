using Aqua
using ElectromagneticFields
using Test

# Package-level quality assurance: type piracy, method ambiguities, stale and duplicated
# dependencies, undefined exports, unbound type parameters, `Project.toml` validity. These are the
# faults the rest of the suite is structurally unable to see — it exercises fields and the values
# they generate, and every one of these is a property of the package as a whole.
#
# Three of them are live concerns here.
#
# `piracies` covers the eight methods this package adds to `GeometricBase`'s generics: `periodic` on
# five equilibrium types, and `periodic`, `functions` and `parameters` on `FieldFunctions`. Each one
# dispatches on a type this package owns, which is extension rather than piracy, and this is what
# says so. A new chart family brings a `periodic` method with it, so the property is worth asserting
# where a list of the sites would go stale.
#
# `undefined_exports` covers an export list dominated by the musical-isomorphism accessors, which
# are defined across the files of `src/analytic/` rather than in the module that exports them. A
# rename in one place and not the other is invisible until a caller reaches for the name.
#
# `stale_deps` covers `ConstructionBase`, which `src/` reaches only through qualified access.
# `fatou lint` reads that as an unused import; Aqua resolves the qualified call and does not.
Aqua.test_all(ElectromagneticFields)
