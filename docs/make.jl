using NOVAS
using Documenter

DocMeta.setdocmeta!(NOVAS, :DocTestSetup, :(using NOVAS); recursive=true)

makedocs(;
    modules=[NOVAS],
    authors="Kiran Shila <me@kiranshila.com> and contributors",
    repo="https://github.com/kiranshila/NOVAS.jl/blob/{commit}{path}#{line}",
    sitename="NOVAS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kiranshila.github.io/NOVAS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Benchmarks" => "benchmarks.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/kiranshila/NOVAS.jl",
    devbranch="main",
)
