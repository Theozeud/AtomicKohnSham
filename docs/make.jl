push!(LOAD_PATH, "../src/")

using Documenter, AtomicKohnSham

makedocs(sitename = "AtomicKohnSham Documentation",
    pages = [
        "index.md",
        "Tutorials" => [
            "Models" => "tutorials/model.jl",
            "AtomProblem" => "tutorials/problem.jl",
            "Groundstate" => "tutorials/groundstate.jl",
            "Analysing results" => "tutorials/solution.jl"
        ]
    ])

deploydocs(
    repo = "github.com/Theozeud/AtomicKohnSham.git",
)
