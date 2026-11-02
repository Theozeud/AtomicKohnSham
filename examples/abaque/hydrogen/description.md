# Description

First in the `examples/abaque` series (`oxygen`, `Z20`, `Z50` follow the
same code shape, generalized to heavier atoms -- see their own
`description.md` for what changes and why).

An abaque: `|Etot - Etot_ref|` vs `Rmax` (log-log), for hydrogen (`Z=N=1`),
`nspin=1`, **Slater exchange only** (`BuiltinFunctional(:lda_x; nspin=1)`,
`ec = NoFunctional(1)` -- no correlation). One curve per `(Nmesh, ordermax)`
combination: `Nmesh in [5, 10, 20]` via color (a single-hue sequential ramp,
light-to-dark blue, since `Nmesh` is an ordered quantity, not an unordered
category), `ordermax in [2, 5, 10]` via line style (solid/dash/dot).

`lh`, `nh` sized via the same `basis_size(Z)` heuristic used throughout the
other examples (`Z=1` -> `lh=1, nh=2`).

## Why `ordermax` stops at 10

`ordermax=20` was tried first and dropped: `P1IntLegendreBasis` hits
ill-conditioned mass matrices at high polynomial degree -- confirmed
directly, `ordermax=24` throws `InexactError` from complex eigenvalues even
at a small `Nmesh=10`, and `ordermax=22` stops converging cleanly at all
(`MAXITERS` after 400 iterations, 193s, at `Nmesh=50`). This matches the
open "Why is there a problem when taking polynomials with too big degree?"
item in the top-level `ToDo` file -- a real, pre-existing limitation of this
basis, not something specific to this script. `ordermax=10` is comfortably
inside the stable range and already near-exact for hydrogen on its own:
`examples/discretization convergence` shows `ordermax=9` alone reaching
`~1e-14`.

## Reference

`Etot_ref` comes from one reference run at `Rmax=2000, Nmesh=80,
ordermax=13, scftol=1e-11` -- a clear notch beyond every grid corner in all
three knobs at once. `run.jl` checks this automatically against
`target_err = 1e-10`: the best grid corner reaches `2.9e-15`, so the
reference is comfortably resolved.

Grid solves use `scftol=1e-10` -- tighter values (`1e-11`/`1e-12`) were
tried first and made the SCF stall (`MAXITERS`) at larger `Nmesh`/`ordermax`
corners, apparently against a floating-point noise floor that rises with
degrees of freedom, not a real physical limit.

## Results (`abaque.txt`, `abaque_hydrogen.pdf`)

`abaque.txt` has one `(Nmesh, Rmax)` table block per `ordermax`. The plot
overlays all nine curves at once and makes the interaction between the two
knobs immediately visible:

- **`ordermax=2` (solid) plateaus around `1e-2`-`1e0`, almost independent of
  `Nmesh`** -- at this degree the FEM basis simply can't represent the
  orbital's curvature, no matter how many elements it's cut into. More
  `Nmesh` barely moves these three lines.
- **`ordermax=5` (dashed) drops to `~1e-6`-`1e-9`**, with real separation
  between `Nmesh` values starting to show.
- **`ordermax=10` (dotted) reaches `~1e-14`-`1e-15`**, and this is where
  `Nmesh` matters most: `Nmesh=20` clears `Nmesh=5` by several orders of
  magnitude once `Rmax` is large enough.
- **Every curve still collapses together for `Rmax <= 10`**, same as the
  2-knob version of this study -- domain truncation dominates everyone's
  error there regardless of `Nmesh` or `ordermax`.

The takeaway: `ordermax` sets the *ceiling* on achievable accuracy (you
cannot out-refine-`Nmesh` your way past a too-low polynomial degree), while
`Nmesh` is what lets you approach that ceiling. This is the p-refinement vs
h-refinement relationship from `examples/discretization convergence`, now
visible directly rather than inferred from two separate 1D sweeps.

## Running it

`run.jl` activates the shared `examples/` environment itself, so just run:
```julia
julia run.jl
```
72 grid points (`3 Nmesh x 3 ordermax x 8 Rmax`) + 1 reference, all with a
small basis (hydrogen, `lh=1, nh=2`) -- takes about a minute.
