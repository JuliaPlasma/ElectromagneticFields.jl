module ElectromagneticFields

using GeometricBase
using LinearAlgebra
using PrecompileTools
using RuntimeGeneratedFunctions
using StaticArrays
using Symbolics

import ConstructionBase
import NaNMath

import GeometricBase: functions, parameters, periodicity

# This module is the default `cache_module` for the code `FieldFunctions` generates, so the bodies
# built during the precompile workload below land here and survive into the package image.
RuntimeGeneratedFunctions.init(@__MODULE__)

export ElectromagneticField

include("field.jl")

export AnalyticField, AnalyticEquilibrium, AnalyticPerturbation, ZeroPerturbation
export CartesianField, CartesianEquilibrium, CartesianPerturbation

# The ITER parameters shared by the three axisymmetric tokamak charts, which describe the same
# device in different coordinates.
const ITER_R₀ = 6.2
const ITER_B₀ = 5.3
const ITER_q₀ = √2

include("analytic/analytic_field.jl")
include("analytic/cartesian_field.jl")

export FieldFunction, FieldFunctions, @precompilable_fields, clear_field_cache!
export functions, parameters, periodicity, coordinates, orientation
export equilibrium, perturbation
export to_cartesian, from_cartesian, DF, DF̄, J, rangemin, rangemax
export g♭, g♯, Dg♭, Dg♯, DDg♭, DDg♯
export A♭, A♯, DA♭, DDA♭, φ
export B, DB, DDB, B♭, B♯, B♮, B♭♭, DB♭
export b♭, b♯, b♮, Db♭, Db♮, DDb♭
export a♭, a♯, a♮, c♭, c♯, c♮
export E♭, E♯, DE♭

include("analytic/field_functions.jl")

export ABCEquilibrium
export AxisymmetricTokamakCartesianEquilibrium, AxisymmetricTokamakCartesianITER
export AxisymmetricTokamakCylindricalEquilibrium, AxisymmetricTokamakCylindricalITER
export AxisymmetricTokamakToroidalEquilibrium, AxisymmetricTokamakToroidalITER
export AxisymmetricTokamakToroidalRegularizationEquilibrium
export DipoleField
export EzCosZPerturbation
export PenningTrapAsymmetricEquilibrium
export PenningTrapBottleEquilibrium
export PenningTrapUniformEquilibrium
export QuadraticPotentialsField
export AbstractSolovevEquilibrium
export SolovevEquilibrium, SolovevEquilibriumITER, SolovevEquilibriumNSTX,
       SolovevEquilibriumFRC
export SolovevXpointEquilibrium, SolovevXpointEquilibriumITER, SolovevXpointEquilibriumNSTX
export SolovevDoubleXpointEquilibrium, SolovevDoubleXpointEquilibriumNSTX
export SolovevSymmetricEquilibrium
export SingularEquilibrium
export SymmetricQuadraticEquilibrium
export ThetaPinchEquilibrium

include("analytic/abc.jl")
include("analytic/axisymmetric_tokamak_cartesian.jl")
include("analytic/axisymmetric_tokamak_cylindrical.jl")
include("analytic/axisymmetric_tokamak_toroidal.jl")
include("analytic/axisymmetric_tokamak_toroidal_regularization.jl")
include("analytic/dipole.jl")
include("analytic/ezcosz.jl")
include("analytic/penning_trap_asymmetric.jl")
include("analytic/penning_trap_bottle.jl")
include("analytic/penning_trap_uniform.jl")
include("analytic/quadratic_potentials.jl")
include("analytic/solovev_abstract.jl")
include("analytic/solovev.jl")
include("analytic/solovev_symmetric.jl")
include("analytic/singular.jl")
include("analytic/symmetric_quadratic.jl")
include("analytic/theta_pinch.jl")

export plot_equilibrium, plot_equilibrium!

include("plots.jl")

# Building a field in a fresh session is dominated by Julia compiling Symbolics' own machinery for
# the expression types a trace produces, not by this package: on the toroidal tokamak the symbolic
# trace and code generation account for 11 of the 12 seconds and compiling the generated code for
# 0.6. The cost is per expression shape, and because the generated code takes the parameters as an
# argument rather than baking them in, the shape is fixed by the equilibrium's type alone.
#
# Tracing one field of each type here therefore settles it for every parameter value of that type.
# The specializations land in the package image, and so does `FIELD_CACHE`, so a fresh session
# finds the generated functions already built: all twenty shipped equilibria together cost 0.01 s
# to construct.
#
# The price is this package's own precompilation, about 25 s, paid once per version. A downstream
# package can do the same for a field of its own; see `@precompilable_fields`.
@setup_workload begin
    equilibria = (ABCEquilibrium(),
        AxisymmetricTokamakCartesianEquilibrium(),
        AxisymmetricTokamakCylindricalEquilibrium(),
        AxisymmetricTokamakToroidalEquilibrium(),
        AxisymmetricTokamakToroidalRegularizationEquilibrium(),
        DipoleField(),
        PenningTrapUniformEquilibrium(),
        PenningTrapBottleEquilibrium(),
        PenningTrapAsymmetricEquilibrium(),
        QuadraticPotentialsField(),
        SingularEquilibrium(),
        SolovevEquilibriumITER(),
        SolovevXpointEquilibriumITER(),
        SolovevSymmetricEquilibrium(),
        SymmetricQuadraticEquilibrium(),
        ThetaPinchEquilibrium())
    ξ = [0.5, 0.5, 0.5]

    @compile_workload begin
        for equ in equilibria
            field = FieldFunctions(equ)
            for name in FIELD_FUNCTION_NAMES
                getfield(@__MODULE__, name)(field, 0.0, ξ)
            end
            parameters(field)
            coordinates(field)
            orientation(field)
        end
    end
end

end
