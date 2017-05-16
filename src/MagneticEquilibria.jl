module MagneticEquilibria

    export MagneticEquilibrium

    include("equilibrium.jl")

    export AnalyticEquilibrium, load_equilibrium

    include("analytic/analytic_equilibrium.jl")

    export ABC, AxisymmetricTokamak, Solovev, SolovevXpoint, SymmetricQuadratic
    export analyticA₁, analyticA₂, analyticA₃, analyticMetric

    include("analytic/abc.jl")
    include("analytic/axisymmetric_tokamak.jl")
    include("analytic/solovev_psi.jl")
    include("analytic/solovev.jl")
    include("analytic/solovev_xpoint.jl")
    include("analytic/symmetric_quadratic.jl")

end
