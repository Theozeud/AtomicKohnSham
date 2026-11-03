# Description

ODA against QCFON on neutral atoms (C, Ne, Na, Cl) in LSDA (Slater + PW92), with
the same discretization and tolerance (`scftol = 1e-11`).

Both algorithms are measured with the same quantity, computed at every iteration
independently of the algorithm by a callback: the self-consistency residual
`‖D_out[D] - D‖`, where `D` is the current density matrix and `D_out[D]` is built
from the eigenvectors of `H[D]` with the current occupations. For ODA, `D` is the mixed density
matrix of the line search; for QCFON, it is the density matrix that QCFON keeps
as its current solution, `D_out[D_n]` built from the Newton iterate `D_n` (this
is why the first Newton point can lie above the last ODA one: the output of a
not yet converged density matrix is less converged than the ODA mixture). The
energy error is measured against the QCFON energy.

`run.jl` produces:

- `convergence_iterations.pdf`: residual and energy error vs iteration (dots:
  Newton iterations of QCFON, after the ODA warm-up);
- `convergence_time.pdf`: residual vs wall time, excluding compilation and the
  time spent in the monitoring callback (fastest of 3 runs);
- `comparison.txt`: iterations, times, final residuals and energy differences.

`run.jl` activates the shared `examples/` environment itself, so just run:
```julia
julia run.jl
```
