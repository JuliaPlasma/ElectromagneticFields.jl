using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Plots" include("plots.jl")
end
if "slow" in GROUPS
    @safetestset "Analytic equilibria" include("analytic/equilibria.jl")
end
