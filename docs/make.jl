using Documenter, ElectromagneticFields

makedocs(
    sitename = "ElectromagneticFields.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["Home" => "index.md",
             "Modules" => "modules.md",
            ]
)

deploydocs(
    repo   = "github.com/DDMGNI/ElectromagneticFields.jl.git"
)
