using Documenter, ElectromagneticFields
using CairoMakie

CairoMakie.activate!(type="png")

makedocs(
    sitename="ElectromagneticFields.jl",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        assets=[asset("assets/style.css", class=:css, islocal=true)]),
    pages=["Home" => "index.md",
        "Usage" => "usage.md",
        "Coordinates" => "coordinates.md",
        "Fields" => "fields.md",
        "Plotting" => "plotting.md",
        "Analytic Fields" => [
            "Overview" => "analytic/index.md",
            "Arnold-Beltrami-Childress Field" => "analytic/abc.md",
            "Axisymmetric Tokamak (Cartesian)" => "analytic/axisymmetric_tokamak_cartesian.md",
            "Axisymmetric Tokamak (Cylindrical)" => "analytic/axisymmetric_tokamak_cylindrical.md",
            "Axisymmetric Tokamak (Toroidal)" => "analytic/axisymmetric_tokamak_toroidal.md",
            "Dipole" => "analytic/dipole.md",
            "Penning Traps" => "analytic/penning_traps.md",
            "Quadratic Potentials" => "analytic/quadratic_potentials.md",
            "Singular Field" => "analytic/singular.md",
            "Solov'ev Equilibrium" => "analytic/solovev.md",
            "Symmetric Solov'ev Equilibrium" => "analytic/solovev_symmetric.md",
            "Symmetric Quadratic Field" => "analytic/symmetric_quadratic.md",
            "Theta Pinch" => "analytic/theta_pinch.md",
        ],
        "Modules" => "modules.md",
    ]
)

deploydocs(
    repo="github.com/JuliaPlasma/ElectromagneticFields.jl"
)
