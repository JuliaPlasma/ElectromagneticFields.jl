
using ElectromagneticFields
using SafeTestsets
using Test

@safetestset "Aqua" begin
    include("aqua_tests.jl")
end

include("test_analytic.jl")
include("test_plots.jl")
