module ElectromagneticFields

    using Documenter
    using LinearAlgebra

    export MagneticEquilibrium

    include("equilibrium.jl")

    export AnalyticEquilibrium, AnalyticPerturbation, ZeroPerturbation
    export load_equilibrium, periodicity

    include("analytic/analytic_equilibrium.jl")

    export ABC, AxisymmetricTokamakCartesian, AxisymmetricTokamakCylindrical,
           AxisymmetricTokamakToroidal, EzCosZ,
           Solovev, SolovevXpoint, SolovevQuadratic, SymmetricQuadratic, ThetaPinch
    export SolovevITER, SolovevNSTX, SolovevXpointITER, SolovevXpointNSTX
    export @abc_equilibrium,
           @axisymmetric_tokamak_equilibrium_cartesian,
           @axisymmetric_tokamak_equilibrium_cylindrical,
           @axisymmetric_tokamak_equilibrium_toroidal,
           @ezcosz_perturbation,
           @solovev_equilibrium,
           @solovev_xpoint_equilibrium,
           @solovev_equilibrium_quadratic,
           @symmetric_quadratic_equilibrium,
           @theta_pinch_equilibrium

    include("analytic/abc.jl")
    include("analytic/axisymmetric_tokamak_cartesian.jl")
    include("analytic/axisymmetric_tokamak_cylindrical.jl")
    include("analytic/axisymmetric_tokamak_toroidal.jl")
    include("analytic/ezcosz.jl")
    include("analytic/solovev_common.jl")
    include("analytic/solovev_psi.jl")
    include("analytic/solovev.jl")
    include("analytic/solovev_quadratic.jl")
    include("analytic/solovev_xpoint.jl")
    include("analytic/symmetric_quadratic.jl")
    include("analytic/theta_pinch.jl")

end
