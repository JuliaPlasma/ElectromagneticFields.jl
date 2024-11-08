using Documenter, ElectromagneticFields

ENV["GKSwstype"] = "100"

makedocs(
    sitename = "ElectromagneticFields.jl",
    format = Documenter.HTML(
                prettyurls = get(ENV, "CI", nothing) == "true",
                assets = [asset("assets/style.css", class=:css, islocal=true)]),
    pages = ["Home" => "index.md",
             "Analytic Fields" => "analytic.md",
             "Modules" => "modules.md",
            ]
)

deploydocs(
    repo   = "github.com/JuliaPlasma/ElectromagneticFields.jl"
)
