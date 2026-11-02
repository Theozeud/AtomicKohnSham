# Description

Fourth in the `examples/abaque` series (see `hydrogen`, `oxygen`, `Z20`).
Abaque: `|Etot - Etot_ref|` vs `Rmax` (log-log), for `Z=50` (tin,
`Z=N=50`), `nspin=1`, **Slater exchange only** (`BuiltinFunctional(:lda_x;
nspin=1)`, `ec = NoFunctional(1)` -- no correlation). Same grid and code
shape as oxygen/`Z20`: `Nmesh in [8, 16, 32]` via color, `ordermax in
[2, 5, 10]` via line style.

`lh`, `nh` sized via the same `basis_size(Z)` heuristic (`Z=50` ->
`lh=3, nh=10`).

## Two things that needed adjusting for tin specifically

- **`Rmax_ref = 1500`, not `2000`.** At `Rmax=2000` the SCF stalls at
  `MAXITERS` for tin even with a generous, fixed `(Nmesh, ordermax)` --
  checked directly: 300 iterations, still not converged. The same
  discretization at `Rmax=1500` converges cleanly in under 100 iterations.
  Not chased further -- `1500` already sits comfortably beyond the grid's
  own `Rmax=1000` ceiling, which is all the reference needs to do.
- **`target_err = 1e-2`**, the loosest in the series. Like calcium (`Z20`),
  the grid's own `Rmax=1000` ceiling carries real truncation error here --
  `~3.5e-3` Ha for tin, roughly an order of magnitude worse than calcium's
  `~2.8e-5`. The trend across the series (`hydrogen -> oxygen -> Z20 ->
  Z50`) is monotonic: heavier atoms need a bigger `Rmax` to reach the same
  relative accuracy, and this grid's fixed `Rmax=1000` ceiling increasingly
  becomes the binding constraint rather than `Nmesh`/`ordermax`.

## Reference

`Etot_ref` comes from one reference run at `Rmax=1500, Nmesh=64,
ordermax=16, scftol=1e-8`. Best grid error: `2.0e-3` (comfortably under the
`1e-2` target -- consistent with the `Rmax=1000` ceiling estimate above).

## Running it

```julia
julia run.jl
```
72 grid points + 1 reference -- the slowest in the series (tin's SCF needs
more iterations to converge at several grid points), but still well under
a couple of minutes.
