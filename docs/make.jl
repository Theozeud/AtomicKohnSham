push!(LOAD_PATH, "../src/")

using Documenter, AtomicKohnSham

makedocs(sitename = "AtomicKohnSham.jl",
    authors = "Théo Duez",
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "The physical model" => "tutorials/model.md",
            "Discretization" => "tutorials/discretization.md",
            "Solving for the ground state" => "tutorials/groundstate.md",
            "Analyzing a solution" => "tutorials/solution.md",
            "Reports, logs & plots" => "tutorials/results.md",
        ],
        "API Reference" => "reference.md",
    ])

deploydocs(
    repo = "github.com/Theozeud/AtomicKohnSham.git",
)
