module MagneticEquilibria

    export MagneticEquilibrium

    include("equilibrium.jl")

    export AnalyticEquilibrium, load_equilibrium

    include("analytic/analytic_equilibrium.jl")

    export ABC, AxisymmetricTokamak, AxisymmetricTokamakToroidal,
           Solovev, SolovevXpoint, SymmetricQuadratic
    export analyticA₁, analyticA₂, analyticA₃, analyticMetric
    export @axisymmetric_tokamak_equilibrium,
           @axisymmetric_tokamak_equilibrium_toroidal,
           @solovev_equilibrium,
           @solovev_xpoint_equilibrium,
           @symmetric_quadratic_equilibrium

    include("analytic/abc.jl")
    include("analytic/axisymmetric_tokamak_cartesian.jl")
    include("analytic/axisymmetric_tokamak_toroidal.jl")
    include("analytic/solovev_psi.jl")
    include("analytic/solovev.jl")
    include("analytic/solovev_xpoint.jl")
    include("analytic/symmetric_quadratic.jl")

end
