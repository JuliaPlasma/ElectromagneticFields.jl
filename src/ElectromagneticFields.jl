module ElectromagneticFields

using GeometricBase
using LinearAlgebra
using StaticArrays
using Symbolics

import NaNMath

import GeometricBase: functions, parameters, periodicity

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

export FieldFunction, FieldFunctions
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

end
