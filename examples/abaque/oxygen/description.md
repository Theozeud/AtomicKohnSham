# Description

Second in the `examples/abaque` series (see `examples/abaque/hydrogen` for
the first). Abaque: `|Etot - Etot_ref|` vs `Rmax` (log-log), for oxygen
(`Z=N=8`), `nspin=1`, **Slater exchange only** (`BuiltinFunctional(:lda_x;
nspin=1)`, `ec = NoFunctional(1)` -- no correlation). One curve per
`(Nmesh, ordermax)` combination: `Nmesh in [8, 16, 32]` via color
(single-hue sequential ramp, light-to-dark blue), `ordermax in [2, 5, 10]`
via line style (solid/dash/dot).

`lh`, `nh` sized via the same `basis_size(Z)` heuristic used throughout
(`Z=8` -> `lh=1, nh=8`).

## Two differences from the hydrogen version

- **`Nmesh_values` starts at 8, not 5.** `basis_size(8)` needs `nh=8`
  orbitals per l-channel, and hydrogen's smallest grid point
  (`Nmesh=5, ordermax=2`, `dof=7`) is too small to even hold that many --
  it throws `ArgumentError("illegal choice of eigenvalue indices...")`
  rather than a merely-inaccurate result. `Nmesh=8, ordermax=2` (`dof=13`)
  clears the floor with room to spare.
- **`target_err = 1e-8`, not `1e-10`.** An 8-electron SCF doesn't settle as
  tightly as hydrogen's 1-electron one at a fixed `scftol`. The best grid
  corner reaches `6.9e-10` here regardless -- comfortably under target, and
  actually not far off hydrogen's own precision.

`ordermax` still stops at 10, same reason as hydrogen (`ordermax=20+` is
numerically unstable for `P1IntLegendreBasis` -- see hydrogen's `run.jl`).

## Reference

`Etot_ref` comes from one reference run at `Rmax=2000, Nmesh=48,
ordermax=13, scftol=1e-9` -- comfortably beyond every grid corner. Best grid
error: `6.9e-10`.

## Results

Same shape as hydrogen's abaque, shifted up in absolute error (more
electrons, more self-interaction to resolve) but with the same structure:
`ordermax=2` plateaus almost independent of `Nmesh`, `ordermax=10` reaches
the lowest floor, and all nine curves collapse together for small `Rmax`
where domain truncation dominates regardless of the other two knobs.

## Running it

```julia
julia run.jl
```
72 grid points + 1 reference -- takes about a minute.
