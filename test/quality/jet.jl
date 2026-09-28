using ElectromagneticFields
using JET
using Test

using ElectromagneticFields: FIELD_FUNCTION_NAMES

# JET's optimisation analysis of the hot path: the generated field functions, which
# `analytic/equilibria.jl` asserts allocation free. Each accessor is analysed at the concrete
# argument types that file passes — a `FieldFunctions` of each equilibrium, a `Float64` time and a
# `Vector{Float64}` of coordinates. A runtime dispatch on that path is what would make a call
# allocate, and this reports it at the line that causes it.
#
# JET does not support every Julia version (the `pre` and `nightly` jobs); there it records one
# `@test_skip`.

const JET_WORKS = isdefined(JET, :JET_AVAILABLE) ? JET.JET_AVAILABLE : JET.JET_LOADABLE

# The constructor names label the testsets: several constructors return the same type.
const EQUILIBRIA = (
    ("ABCEquilibrium", ABCEquilibrium()),
    ("AxisymmetricTokamakCartesianEquilibrium", AxisymmetricTokamakCartesianEquilibrium()),
    ("AxisymmetricTokamakCylindricalEquilibrium",
        AxisymmetricTokamakCylindricalEquilibrium()),
    ("AxisymmetricTokamakToroidalEquilibrium", AxisymmetricTokamakToroidalEquilibrium()),
    ("AxisymmetricTokamakToroidalRegularizationEquilibrium",
        AxisymmetricTokamakToroidalRegularizationEquilibrium()),
    ("DipoleField", DipoleField()),
    ("PenningTrapUniformEquilibrium", PenningTrapUniformEquilibrium()),
    ("PenningTrapBottleEquilibrium", PenningTrapBottleEquilibrium()),
    ("PenningTrapAsymmetricEquilibrium", PenningTrapAsymmetricEquilibrium()),
    ("QuadraticPotentialsField", QuadraticPotentialsField()),
    ("SingularEquilibrium", SingularEquilibrium()),
    ("SymmetricQuadraticEquilibrium", SymmetricQuadraticEquilibrium()),
    ("ThetaPinchEquilibrium", ThetaPinchEquilibrium()),
    ("SolovevEquilibriumFRC", SolovevEquilibriumFRC()),
    ("SolovevEquilibriumITER", SolovevEquilibriumITER()),
    ("SolovevXpointEquilibriumITER", SolovevXpointEquilibriumITER()),
    ("SolovevEquilibriumNSTX", SolovevEquilibriumNSTX()),
    ("SolovevXpointEquilibriumNSTX", SolovevXpointEquilibriumNSTX()),
    ("SolovevDoubleXpointEquilibriumNSTX", SolovevDoubleXpointEquilibriumNSTX()),
    ("SolovevSymmetricEquilibrium", SolovevSymmetricEquilibrium())
)

if JET_WORKS
    @testset "$label" for (label, equ) in EQUILIBRIA
        field = FieldFunctions(equ)
        @testset "$name" for name in FIELD_FUNCTION_NAMES
            f = getfield(ElectromagneticFields, name)
            @test isempty(JET.get_reports(JET.report_opt(f,
                (typeof(field), Float64, Vector{Float64});
                target_modules = (ElectromagneticFields,))))
        end
    end
else
    @test_skip "JET does not work on this Julia version"  # aviatesk/JET.jl#681
end
