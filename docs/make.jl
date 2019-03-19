using Documenter, ElectromagneticFields

makedocs(
    sitename = "ElectromagneticFields.jl",
    format = Documenter.HTML(),
    pages = ["Home" => "index.md",
             "Modules"    => "modules.md",
             ]
)

deploydocs(
    repo   = "github.com/DDMGNI/ElectromagneticFields.jl.git"
)
