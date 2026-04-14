using Documenter
using SynthDiD

makedocs(
    sitename = "SynthDiD.jl",
    modules = [SynthDiD],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Introduction"      => "tutorials/01_introduction.md",
            "More on Plotting"  => "tutorials/02_more_plotting.md",
            "Paper Results"     => "tutorials/03_paper_results.md",
        ],
        "Reference" => "reference.md",
    ],
    warnonly = true,
)
