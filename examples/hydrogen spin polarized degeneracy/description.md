# Description

Hydrogen, `Z=N=1`, spin-polarized (`nspin=2`), `lh=0`/`nh=1` (1s only) --
the smallest possible system that isolates a single twofold degeneracy: the
`1sUP`/`1sDOWN` pair competing for exactly 1 electron. Built to test the
`flat_tol` guard added to `resolve_twofold_degeneracy!`
(`src/algorithms/aufbau/optimized.jl`), which falls back to a
degeneracy-weighted symmetric occupation split whenever the total energy
sampled across the occupation-mixing line is flat to within `flat_tol`,
instead of trusting `line_search_energy`'s minimizer (which is exact when
there's real curvature to follow, but becomes an artifact of roundoff once
that curvature collapses to zero).

Same discretization for every case (`Rmax=1000`, `Nmesh=30`, `ordermax=10`,
`ODA(tinit=0.6, scftol=1e-10)`, `maxiter=300`). Five cases:

| case                          | functional              | aufbau                              |
|--------------------------------|---------------------------|----------------------------------------|
| `hartree_only_guard_on`         | none (Hartree only)         | `OptimizedAufbau(max_degen=2)` (default `flat_tol=1e-8`) |
| `hartree_only_guard_off`        | none (Hartree only)         | `OptimizedAufbau(max_degen=2, flat_tol=0)` -- guard disabled, old behavior |
| `lda_guard_on`                  | Slater + Perdew-Wang (Libxc) | `OptimizedAufbau(max_degen=2)`          |
| `lda_guard_off`                 | Slater + Perdew-Wang (Libxc) | `OptimizedAufbau(max_degen=2, flat_tol=0)` |
| `hartree_only_smeared`          | none (Hartree only)         | `SmearedAufbau(Temp=0.01)` -- the alternative, smoothing-based route |

Each case logs `log.txt`/`results.txt`, plus a per-iteration occupation
tracker (`OccTracker`, a small custom `callback!` defined right in
`run.jl`) records `n(1sUP)` at every SCF step so the trajectory -- not just
the final answer -- can be plotted (`occupation_comparison.pdf`).

## What was actually observed

```
SCF OUTCOME
Hartree only, guard ON       SUCCESS in 16 iterations
Hartree only, guard OFF      SUCCESS in 16 iterations
LDA, guard ON                SUCCESS in 14 iterations
LDA, guard OFF               SUCCESS in 14 iterations
Hartree only, SmearedAufbau  SUCCESS in 16 iterations

FINAL OCCUPATION (n(1sUP), n(1sDOWN))
Hartree only, guard ON       n_up = 0.5   n_down = 0.5
Hartree only, guard OFF      n_up = 0.0   n_down = 1.0
LDA, guard ON                n_up = 0.0   n_down = 1.0
LDA, guard OFF               n_up = 0.0   n_down = 1.0
Hartree only, SmearedAufbau  n_up = 0.5   n_down = 0.5
```

This is **not** what the earlier discussion predicted, and it's worth being
upfront about that rather than papering over it: `occupation_comparison.pdf`
shows every case converging smoothly and monotonically to a fixed
occupation by iteration 2, and `convergence_comparison.pdf` shows every
stopping-criterion curve (including both `guard OFF` runs) decaying
cleanly to `scftol` -- no oscillation, no stalling, `guard OFF` reaches
`SUCCESS` in exactly the same iteration count as `guard ON`.

**Why the predicted oscillation didn't show up here.** The oscillation
mechanism described earlier relies on `ε(1sUP)` and `ε(1sDOWN)` being
computed independently and only *approximately* equal, so that roundoff
noise in that near-tie flips sign unpredictably between iterations. For
`hartree_only`, that's not what happens: the Hartree potential is built
once from the *total* density and reused bit-for-bit for both spin
channels' Hamiltonians (`assemble.jl`: `vH = Hfix[l+1] + VxcUP + Hartree`
and `... + VxcDOWN + Hartree`, same `Hartree` array both times, and
`VxcUP=VxcDOWN=0` identically since there's no XC functional here at all).
`find_orbital!` therefore diagonalizes the *exact same matrix* for both spin
channels, every iteration -- `ε(1sUP)` and `ε(1sDOWN)` are equal to the
last bit, not just numerically close, and so are the two extreme densities
`D1`/`D2` built from them. The curvature `a` and slope `b` in
`line_search_energy` are exactly `0.0`, not a roundoff-sized residual --
there is no noise here to drive an iteration-to-iteration sign flip.

**What actually goes wrong without the guard, then:** with `a == 0.0`
exactly, the *old* code's tie-break (`line_search_energy`'s
`nl=false` branch: `f0 <= f1 ? 0 : 1`) resolves to `t=0` **deterministically,
every time** -- not randomly. It always dumps the electron entirely into
whichever of the two indices happens to sort second (`1sDOWN`, an artifact
of `collect_degenerate_block!`'s scan order over `indices_sort`, itself an
artifact of how `convert_index` lays spin out in the flattened orbital
index -- nothing physical). That's a real bug, just a quieter one than
predicted: it's stable (so it doesn't break convergence), but it silently
reports a **fully spin-polarized state as if it were the answer**, for a
system where nothing in the physics prefers `DOWN` over `UP` -- the
"choice" is coming from array layout, not energetics. `Hartree only, guard
ON` and `Hartree only, SmearedAufbau` -- two independent mechanisms --
agree on the actual least-biased answer (`n=0.5/0.5`), which the old
deterministic tie-break simply never considers.

**The LDA cases show the fix is a no-op when a real minimum exists**, as
intended: `guard ON`/`guard OFF` land on the *identical* answer
(`n_up=0, n_down=1`, `Etot=-0.478711 Ha` both), because Slater exchange
gives a single electron split 50/50 across spin channels strictly higher
exchange energy than the same electron fully in one channel (a real,
`~0.01-0.1 Ha`-scale energy difference here, nowhere near `flat_tol`'s
`1e-12`), so the flatness check correctly declines to override
`line_search_energy` -- which correctly finds the true (if two-fold, by the
`UP<->DOWN` mirror symmetry) minimum.

## Where genuine oscillation *would* show up

The bit-exact tie above is a special property of `hartree_only`'s
spin-blind Hamiltonian construction (and of hydrogen's single electron,
which makes the two `line_search_energy` endpoints literally the same
real-space density, just relabeled). It's the cleanest possible case to
*prove* the fix does the right thing, but it doesn't reproduce the noisier
failure mode motivating it. That mode needs a near-degeneracy between
orbitals that are **not** related by an exact symmetry -- so their energies
really are computed via independent numerical paths and can drift apart or
together by roundoff/discretization-sized amounts iteration to iteration.
That's exactly the `3d`/`4s` near-degeneracy already documented in
`examples/atomic density comparison/run.jl` for iron (`nspin=1`, a
different-`l` tie rather than a spin-mirror one) -- reproducing that same
instability in a spin-polarized, near-but-not-exactly-degenerate open shell
would be the natural next step if it's worth pinning down further.

## Running it

`run.jl` activates the shared `examples/` environment itself
(`Pkg.activate(joinpath(@__DIR__, ".."))`), so just run:
```julia
julia run.jl
```
