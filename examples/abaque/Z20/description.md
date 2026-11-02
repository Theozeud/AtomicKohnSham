# Description

Third in the `examples/abaque` series (see `examples/abaque/hydrogen` and
`examples/abaque/oxygen`). Abaque: `|Etot - Etot_ref|` vs `Rmax` (log-log),
for `Z=20` (calcium, `Z=N=20`), `nspin=1`, **Slater exchange only**
(`BuiltinFunctional(:lda_x; nspin=1)`, `ec = NoFunctional(1)` -- no
correlation). Same grid and code shape as oxygen's: `Nmesh in [8, 16, 32]`
via color, `ordermax in [2, 5, 10]` via line style.

`lh`, `nh` sized via the same `basis_size(Z)` heuristic (`Z=20` ->
`lh=1, nh=10`).

## The one real difference: `Rmax=1000` is not enough for calcium

This is the finding that shaped this folder's parameters. Held at a
generous, fixed `(Nmesh, ordermax)`, `Etot` at `Rmax=1000` differs from
`Rmax=2000` by `~2.8e-5` Ha -- for hydrogen, the equivalent gap was already
down at `~1e-14` by `Rmax=50`. Domain truncation genuinely matters more for
calcium than for the lighter atoms in this series, consistent with a more
diffuse valence density.

Consequence: `target_err = 1e-4`, not `1e-8`/`1e-10` like the lighter
atoms. No `(Nmesh, ordermax)` on this grid can get the error below the
grid's own `Rmax=1000` ceiling (`~2.2e-5` here), no matter how good the
reference is -- that's visible directly in `abaque_Z20.pdf` as the
`ordermax=10` curves flattening out past `Rmax~20` instead of continuing to
drop. `run.jl`'s STEP-4 check is worded to allow for this: a `best_err`
above target can mean either "reference needs pushing further" or "the
grid's own `Rmax` ceiling is the real limit," and this folder is an example
of the second case, confirmed by comparing against `Rmax_ref=2000`.

## Reference

`Etot_ref` comes from one reference run at `Rmax=2000, Nmesh=48,
ordermax=13, scftol=1e-9`. Best grid error: `2.2e-5` (see note above for why
that's the real ceiling, not a reference-quality issue).

## Running it

```julia
julia run.jl
```
72 grid points + 1 reference -- takes about a minute.
