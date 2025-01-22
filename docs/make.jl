using TPSAInterface
using Documenter

DocMeta.setdocmeta!(TPSAInterface, :DocTestSetup, :(using TPSAInterface); recursive=true)

makedocs(;
    modules=[TPSAInterface],
    authors="Matt Signorelli <mgs255@cornell.edu> and contributors",
    sitename="TPSAInterface.jl",
    format=Documenter.HTML(;
        canonical="https://bmad-sim.github.io/TPSAInterface.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bmad-sim/TPSAInterface.jl",
    devbranch="main",
)
