using GeneDrive
using Documenter

makedocs(;
    modules=[GeneDrive],
    authors="Valeri Vasquez",
    repo="https://github.com/vnvasquez/GeneDrive.jl/blob/{commit}{path}#L{line}",
    sitename="GeneDrive.jl",
    format=Documenter.HTML(; prettyurls=prettyurls = get(ENV, "CI", nothing) == "true"),
    pages=["Home" => "index.md", 
    "API Reference" => "api.md", 
    "Examples" => Any["datasetup_examples.md", # Data/getting started 
        "dynamic_examples.md", # ODE
        "decision_examples.md", # optimization
        "developer_examples.md" # Make your own (advanced)
        ]
    ],
)

deploydocs(; repo="github.com/vnvasquez/GeneDrive.jl")
