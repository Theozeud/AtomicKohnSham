# Description

Basic example of a ground-state calculation for an atom/ion (here, spin-unpolarized
sodium with the Slater (`lda_x`) exchange-only functional, no correlation),
demonstrating the full output pipeline: `run.jl` builds the
model/discretization/algorithm, solves it, and produces:

- `log.txt`: parameters + one line per SCF iteration, written incrementally
  through a `LogFileCallback` (survives interruptions/errors, unlike
  `redirect_stdout`).
- `results.txt`: a final report (parameters, energies, occupied orbitals),
  written by `write_report`.
- `density.pdf`, `orbitals.pdf`, `potentials.pdf`, `convergence.pdf`,
  `energy_breakdown.pdf`: publication-style plots produced by
  `AtomicKohnSham`'s `CairoMakie` extension (loaded automatically once
  `using CairoMakie` is added, as done here).

`run.jl` activates the shared `examples/` environment itself
(`Pkg.activate(joinpath(@__DIR__, ".."))`), so just run:
```julia
julia run.jl
```
