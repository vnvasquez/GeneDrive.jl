using GeneDrive
using Documenter

makedocs(;
    modules=[GeneDrive],
    authors="Valeri Vasquez",
    repo="https://github.com/vnvasquez/GeneDrive.jl/blob/{commit}{path}#L{line}",
    sitename="GeneDrive.jl",
    format=Documenter.HTML(; prettyurls=prettyurls = get(ENV, "CI", nothing) == "true", size_threshold = nothing),
    pages=[
        "Home" => "index.md",
        "Features" => "features.md",
        "Tutorials" =>
            Any["datasetup_tutorials.md", "dynamic_tutorials.md", "decision_tutorials.md"],
        #"User Guide" => Any["customization_userguide.md","optimization_userguide.md"]
        "API Reference" => "api.md",
    ],
)

deploydocs(
    repo="github.com/vnvasquez/GeneDrive.jl",
    target="build",
    branch="gh-pages",
    devbranch="main",
    devurl="dev",
    push_preview=true,
    versions=["stable" => "v^", "v#.#"],
)
