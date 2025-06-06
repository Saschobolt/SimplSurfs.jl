using SimplSurfs
using Documenter

DocMeta.setdocmeta!(SimplSurfs, :DocTestSetup, :(using SimplSurfs); recursive=true)

makedocs(;
    modules=[SimplSurfs],
    authors="Sascha Stüttgen <sascha.stuettgen@rwth-aachen.de> and contributors",
    sitename="SimplSurfs.jl",
    format=Documenter.HTML(;
        canonical="https://Saschobolt.github.io/SimplSurfs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Saschobolt/SimplSurfs.jl",
    devbranch="main",
)
