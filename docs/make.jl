push!(LOAD_PATH, "../src/")

using AURORA
using Documenter

makedocs(
        sitename="AURORA.jl",
        modules =[AURORA],
        pages = [
                "Home" => "index.md",
                "Manual" => "manual.md"
                ])

deploydocs(;
    repo="github.com/egavazzi/AURORA.jl",
    )
