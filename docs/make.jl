push!(LOAD_PATH, "../src/")

using Documenter, AtomicKohnSham

makedocs(sitename = "AtomicKohnSham Documentation")

deploydocs(
    repo = "github.com/Theozeud/AtomicKohnSham.git",
)
