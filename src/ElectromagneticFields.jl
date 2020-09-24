__precompile__(false)

module ElectromagneticFields

    using Documenter
    using LinearAlgebra

    export ElectromagneticField

    include("field.jl")

    export AnalyticField, AnalyticEquilibrium, AnalyticPerturbation, ZeroPerturbation

    export load_equilibrium, periodicity

    include("analytic/analytic_equilibrium.jl")
    include("analytic/cartesian_equilibrium.jl")

    using .AnalyticCartesianField

    export ABC, EzCosZ, 
           AxisymmetricTokamakCartesian,
           AxisymmetricTokamakCircular,
           AxisymmetricTokamakCylindrical,
           AxisymmetricTokamakToroidalRegularization,
           Solovev, SolovevXpoint, SolovevQuadratic,
           SymmetricQuadratic, ThetaPinch

    export SolovevITER, SolovevXpointITER,
           SolovevNSTX, SolovevXpointNSTX

    export @abc_equilibrium,
           @axisymmetric_tokamak_equilibrium_cartesian,
           @axisymmetric_tokamak_equilibrium_circular,
           @axisymmetric_tokamak_equilibrium_cylindrical,
           @axisymmetric_tokamak_equilibrium_toroidal_regularisation,
           @ezcosz_perturbation,
           @solovev_equilibrium,
           @solovev_xpoint_equilibrium,
           @solovev_equilibrium_quadratic,
           @symmetric_quadratic_equilibrium,
           @theta_pinch_equilibrium

    include("analytic/abc.jl")
    include("analytic/axisymmetric_tokamak_cartesian.jl")
    include("analytic/axisymmetric_tokamak_circular.jl")
    include("analytic/axisymmetric_tokamak_cylindrical.jl")
    include("analytic/axisymmetric_tokamak_toroidal_regularization.jl")
    include("analytic/ezcosz.jl")
    include("analytic/solovev_abstract.jl")
    include("analytic/solovev.jl")
    include("analytic/solovev_quadratic.jl")
    include("analytic/solovev_xpoint.jl")
    include("analytic/symmetric_quadratic.jl")
    include("analytic/theta_pinch.jl")

end
