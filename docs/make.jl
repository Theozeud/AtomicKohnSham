push!(LOAD_PATH, "../src/")

using Documenter, AtomicKohnSham

makedocs(sitename = "AtomicKohnSham Documentation",
    pages = [
        "index.md",
        "Tutorials" => [
            "Models" => "tutorials/model.md",
            "AtomProblem" => "tutorials/problem.md",
            "Groundstate" => "tutorials/groundstate.md",
            "Analysing results" => "tutorials/solution.md"
        ]
    ])

deploydocs(
    repo = "github.com/Theozeud/AtomicKohnSham.git",
)
