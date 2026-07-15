# Implement GGA (Generalized Gradient Approximation) support

## Context

The code currently only supports LDA/LSDA: the XC potential is a purely
multiplicative operator `v_xc(ρ(r))`, so `assemble_exc!`
(`src/discretization/assemble.jl`) builds it as a single density-weighted
**mass matrix** (`fill_mass_matrix!` with a `FunWeight` closure over `ρ(r)`).
GGA functionals (PBE etc.) additionally depend on `σ = |∇ρ|²`, and their
contribution to the KS potential involves `∇ρ`, not just `ρ` — this needs a
genuinely new FEM building block, not just a different functional plugged
into the existing machinery.

Because this is a **spherically-symmetric radial problem**, GGA collapses to
a much simpler 1D case than general 3D GGA: the density only depends on `r`,
so `∇ρ` is purely radial and `σ = (dρ/dr)²` — no angular cross terms, ever.
This makes the feature tractable without touching the angular/`l`-channel
structure at all.

This pass targets **Libxc GGA functionals only** (PBE via
`Functional(:gga_x_pbe/:gga_c_pbe)`); a from-scratch `BuiltinFunctional` PBE
(mirroring `SlaterXα`/`PerdewWang`) is deferred to a later, separate task.

## The key derivation (do this once, get it right, reuse everywhere)

Let `q(r) = Σᵢⱼ D_ij φᵢ(r)φⱼ(r)` (already computed today as `s` in
`eval_density`/`optimized_eval_density!`), so `ρ(r) = q(r)/(4πr²)`.

- **Density gradient**: `q'(r) = 2 Σᵢⱼ D_ij φᵢ'(r)φⱼ(r)` (using `D`
  symmetric — no need to separately track `φᵢφⱼ'`, halving the work), and
  `ρ'(r) = q'(r)/(4πr²) - 2ρ(r)/r`.
- **XC energy**: unchanged in shape — `Exc = ∫ ρ(r) εxc(ρ,σ) 4πr² dr` via the
  existing quadrature grid, `σ(r) = ρ'(r)²` (nspin=1) or
  `(σ↑↑,σ↑↓,σ↓↓) = (ρup'², ρup'ρdown', ρdown'²)` (nspin=2) fed to Libxc's
  `zk`.
- **XC potential matrix** (the Hamiltonian contribution, `∂Exc/∂D_ij`, via
  the chain rule through `ρ(r)` and `σ(r)`, both linear/quadratic in `D`):

  ```
  ∂Exc/∂D_ij = ∫ [ (vρ - 4vσρ'/r) φᵢφⱼ + 2vσρ' (φᵢ'φⱼ + φᵢφⱼ') ] dr
  ```

  i.e. `V_xc = M[vρ - 4vσρ'/r] + Mmix[2vσρ'] + Mmix[2vσρ']ᵀ`, where `M[w]`
  is the existing weighted mass matrix and `Mmix[w] = ∫w(r)φᵢ'(r)φⱼ(r)dr`
  is a **new** (not-symmetric-by-itself) building block — its sum with its
  own transpose is what's symmetric, matching (and reusing) the existing
  defensive `(V+Vᵀ)/2` step already present in `assemble_exc!`.

  This weak-form/matrix-element route (differentiate `Exc` w.r.t. `D_ij`
  directly) is deliberately used instead of building an explicit real-space
  `v_xc(r) = ∂f/∂ρ - ∇·(2vσ∇ρ)` and integrating it as a mass-matrix weight:
  the latter needs `∇vσ`, which isn't well-defined for the piecewise-`C⁰`
  FEM density. The matrix-element route sidesteps that entirely (standard
  practice in basis-set GGA codes) and reuses 100% of the existing quadrature
  infrastructure.

  Numeric-bookkeeping risk (same class of bug as the PW92 `F2P0`
  factor-of-2 found earlier this session): the factors of 2 and 4 above must
  be re-derived and checked by finite difference before trusting them, not
  just typed in.

## Stage 1 — New FEM primitive: `Mmix[w] = ∫ w(r) φᵢ'(r)φⱼ(r) dr`

Reuses the *exact* existing generator-product / quadrature machinery, just
with a new pairwise product:

- `src/fem/basis.jl` (`BasisCache`): add
  `prodMixedG = mul(derivpolynomials, polynomials)` alongside the existing
  `prodMG = pairwiseproduct(polynomials)` /
  `prodMdG = pairwiseproduct(derivpolynomials)` — `mul` already supports two
  different `PolySet`s (it's how `prodTG = mul(prodMG, polynomials)` is
  built), so no new polynomial-level code is needed.
- `src/fem/integration methods.jl` (`GaussLegendre`): add
  `Qmixedgenx = evaluate(basis.cache.prodMixedG, x)`, mirroring `Qgenx`.
- `src/fem/local matrix.jl`: the
  `fill_local_matrix!(K, ::GaussLegendre, ::FunWeight, eldata, ::PolySet, basis)`
  method currently picks `Qgenx` vs `Pgenx` by **matching `length(K)` against
  their sizes** — this breaks for the mixed case since `Mmix` is also `n×n`
  (same shape as the plain mass matrix). Switch the dispatch to key off
  `eldata.s` (already threaded through, currently `:M`/`:Md`/`:T`) instead of
  shape, adding a new `:Mmix` case that selects `Qmixedgenx`.
- `src/fem/matrices.jl`: add `fill_mixed_mass_matrix!(basis, A; weight, method)`
  mirroring `fill_mass_matrix!`'s per-cell loop (only the "recompute every
  cell" branch is needed — the `NoWeight` reuse-across-cells optimization
  doesn't apply since GGA weights are always density-dependent `FunWeight`s),
  and calls `getelement(basis, k, :Mmix)`.

**Validate** (standalone, before touching anything downstream): for a small
test basis, compare a handful of `Mmix[w]_ij` entries against direct
numerical quadrature (`QuadGK`) of `w(r)φᵢ'(r)φⱼ(r)` for an arbitrary smooth
`w`, the same style of ground-truth check used for `test/fem.jl`'s "FEM
Matrices" testset.

## Stage 2 — Density gradient evaluation

- `src/discretization/density.jl`: add
  `eval_density_gradient!`/`eval_density_gradient`, mirroring
  `eval_density!`/`eval_density` but contracting `D[Ik,Ik]` against the
  mixed-product evaluations instead of `g_ij`, then applying
  `ρ' = q'/(4πr²) - 2ρ/r` (needs `ρ(x)` too, so compute both together in one
  pass to avoid a redundant contraction).
- `src/discretization/assemble.jl`: add `optimized_eval_density_gradient!`,
  mirroring `optimized_eval_density!`'s fast per-cell contraction, but using
  a `Qmixedgenx`-based reshape instead of `Qgenx`.

**Validate**: for a real converged density matrix `D` (e.g. from an existing
Hydrogen/Scandium test), finite-difference `eval_density` at `x±h` and
compare to `eval_density_gradient(x)`.

## Stage 3 — Functional interface & assembly wiring

- `src/physics/models.jl`:
  - Add `is_gga(f) = f.family == :gga` (mirrors existing `is_lda`);
    `has_gga(model) = is_gga(model.exchange) || is_gga(model.correlation)`.
    For this pass, **require both `ex`/`ec` to be the same rung** (both LDA
    or both GGA) — error clearly if mixed, rather than silently doing the
    wrong thing; this covers every practical case (Slater+PW92,
    PBE_X+PBE_C) and defers the mixed-rung edge case indefinitely.
  - Extend `evaluate_zk!`/`evaluate_vrho!` to accept and forward `sigma`,
    and add the `vsigma`-returning variant (needed only when
    `has_gga(model)`); Libxc's own `evaluate!` already dispatches on
    `Val{:gga}` and accepts `sigma=`/`vsigma=` kwargs (confirmed directly
    against the installed `Libxc.jl`), so no changes needed on the Libxc
    side — only in these wrapper functions.
  - Add a `PBE(; Z, N, nspin=1)` convenience constructor (mirrors `Slater`),
    built from `Functional(:gga_x_pbe, n_spin=nspin)` /
    `Functional(:gga_c_pbe, n_spin=nspin)`.
- `src/discretization/assemble.jl` (`assemble_exc!`): when `has_gga(model)`,
  the per-cell weight closures must also compute `ρ'` via
  `optimized_eval_density_gradient!`, build `σ`, call the vσ-returning
  evaluator, and assemble `VxcUP`/`VxcDOWN` as
  `M[vρ - 4vσρ'/r] + Mmix[2vσρ'] + Mmix[2vσρ']ᵀ` (new
  `fill_mixed_mass_matrix!` call added alongside the existing
  `fill_mass_matrix!` call) instead of the pure mass matrix used today.
- `src/discretization/energies.jl` (`compute_exc_energy`): when
  `has_gga(model)`, also evaluate `ρ'` at the global quadrature grid `y`
  (analogous `eval_density_gradient!` call) and pass `sigma` into
  `evaluate_zk!`.

**Validate** (required gate before any SCF run): for a fixed, arbitrary
(not necessarily self-consistent) density matrix `D`, finite-difference
`compute_exc_energy(D + h·e_ij) - compute_exc_energy(D - h·e_ij)` over `2h`
and confirm it matches the assembled `VxcUP[i,j]`/`VxcDOWN[i,j]` entry
directly — this is the standard analytic-vs-numerical-derivative check used
repeatedly this session, and is the most direct way to catch a sign/factor
mistake in the `∂Exc/∂D_ij` formula above before it's laundered through an
SCF loop. Also spot-check that the `-4vσρ'/r` term stays well-behaved (no
blow-up) at the smallest-`r` quadrature nodes on the finest mesh cells.

## Stage 4 — End-to-end test

- Add a GGA testset (`test/physics.jl` or a new `test/gga.jl`), e.g. Hydrogen
  or Helium with `PBE(; Z, N)`, following the existing `run_solve`/
  `check_sanity` pattern from `test/hydrogen.jl`. Assert Tier-1 sanity checks
  (energy balance, electron count, density norm, orthonormality) — **not**
  the virial ratio, since GGA XC (like LDA XC) breaks the pure-Coulomb
  virial identity already documented in `check_sanity`'s docstring.
- Note explicitly (as a known, deliberately deferred gap): `eval_vxc`/
  `eval_effective_potential` in `src/solution/postprocess.jl` and the
  plotting stubs are LDA-only today and won't support GGA post-processing
  until a later pass — not required for the SCF loop itself to work
  correctly, only for diagnostic plotting after the fact.

## Explicitly deferred

- `BuiltinFunctional` PBE (higher-than-Float64 precision GGA) — separate task.
- Mixed-rung models (LDA exchange + GGA correlation or vice versa).
- GGA-aware `eval_vxc`/plotting.

## Status

- Stage 1 (`Mmix` FEM primitive) — done, validated against QuadGK
  (`dev/gga/stage1_mixed_mass_matrix.jl`, all 49 checks pass).
- Stage 2 (density gradient evaluation) — done, validated against finite
  difference and cross-checked general-vs-optimized
  (`dev/gga/stage2_density_gradient.jl`, all checks pass). Found and fixed a
  real bug during validation: `optimized_eval_density_gradient!` was missing
  the physical/reference chain-rule Jacobian factor (`basis.shifts[idxmesh][1]`,
  the same `ϕ[1]` used in `evaluate_deriv!`) on the derivative term — `Qmixedgenx`
  is precomputed at reference-cell quadrature nodes, so a raw contraction is a
  reference-space derivative, not a physical one. `optimized_eval_density!` (no
  derivative) doesn't need this factor, which is why the parallel code for `ρ`
  itself was fine.
- Stage 3 (functional interface & assembly wiring) — done, validated via
  analytic-vs-finite-difference check on both `nspin=1` (PBE/Hydrogen) and
  `nspin=2` (PBE/Lithium) (`dev/gga/stage3_potential_matrix.jl`, all checks
  pass). `is_gga`/`has_gga`, the mixed-rung `KSEModel` constructor check, the
  `PBE(; Z, N, nspin)` constructor, `evaluate_vrho_vsigma!`, and the GGA
  branches of `assemble_exc!`/`compute_exc_energy` are all implemented.
  Found and fixed a second real bug during validation:
  `compute_exc_energy`'s GGA branch initially reused the single global
  Gauss-Legendre grid `y` (fine for LDA, whose ρ-only integrand is C⁰), but
  GGA's `σ=ρ'²` is genuinely discontinuous at every mesh node (this FEM basis
  is only C⁰, so `φᵢ'` jumps across cell boundaries) -- one high-order rule
  spanning all cells converges far too slowly across those jumps to trust an
  analytic-vs-FD check (confirmed by varying `npoints`: the FD ratio crept
  from 0.20 to 0.98 only once `npoints` was pushed ~60x higher). Fixed by
  integrating cell-by-cell with `fem_integration_method`'s per-cell reference
  quadrature, exactly like the matrix assembly already does -- this is now
  the pattern `compute_exc_energy` uses for `has_gga(model)`.
- Stage 4 (end-to-end test) — done: `test/gga.jl` (Hydrogen + Helium PBE
  ground states, Tier-1 sanity checks, mixed-rung rejection), wired into
  `test/runtests.jl`. Dev smoke test at `dev/gga/stage4_endtoend.jl`.
