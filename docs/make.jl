using Documenter
using Graphs
using LeidenClustering

DocMeta.setdocmeta!(LeidenClustering, :DocTestSetup, :(using LeidenClustering); recursive=true)

makedocs(;
    modules=[LeidenClustering],
    sitename="LeidenClustering.jl",
    authors="Eric Feltham",
    format=Documenter.HTML(;
        canonical="https://emfeltham.github.io/LeidenClustering.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => "guide.md",
        "Algorithms" => "algorithms.md",
        "Validation" => "validation.md",
        "API reference" => "api.md",
    ],
    checkdocs=:exports,
    warnonly=false,
)

deploydocs(;
    repo="github.com/emfeltham/LeidenClustering.jl.git",
    devbranch="main",
    push_preview=true,
)
