# Description

One atom, one discretization, five electron-electron interaction models --
to see how much of the total energy, the orbital spectrum, and the density
tail comes from each physical term added on top of the mean-field Hartree
interaction, and to cross-check the pure-Julia `BuiltinFunctional` LDA
against Libxc's reference implementation of the same functionals.

Oxygen, `Z=8, N=8`, spin-unpolarized (`nspin=1`). Same mesh/basis/algorithm
for every case (`Rmax=1000`, `Nmesh=30`, `ordermax=10`, `lh=2`, `nh=5`,
`ODA(tinit=0.6, scftol=1e-10)`, `OptimizedAufbau(max_degen=2)`), so the
comparison isn't confounded by discretization error -- only the
exchange/correlation functional changes between cases (`model.hartree=1`
throughout):

| case             | exchange                     | correlation                   |
|--------------------|--------------------------------|---------------------------------|
| `no_functional`     | `NoFunctional`                  | `NoFunctional`                    |
| `lda_x_only`        | Slater (`BuiltinFunctional`)    | `NoFunctional`                    |
| `lda`               | Slater (`BuiltinFunctional`)    | Perdew-Wang (`BuiltinFunctional`) |
| `lda_libxc`         | `:lda_x` (`Libxc`)              | `:lda_c_pw` (`Libxc`)             |
| `pbe`               | PBE (`Libxc`, GGA)              | PBE (`Libxc`, GGA)                |

Each case gets its own `log.txt`/`results.txt` subfolder (same convention
as `examples/basic_example_doublefloat`). The top-level folder collects the
side-by-side view: `comparison.txt` (energies + occupied-orbital table) and
three overlay plots: `density_comparison.pdf`, `energy_breakdown_comparison.pdf`,
`orbital_energies_comparison.pdf`.

## Results

```
ENERGIES (Ha)
            Hartree only   LDA x-only   LDA          LDA (Libxc)  PBE
Ekin        67.038565      73.925425    74.115171    74.115171    74.734377
Ecou        -166.055729    -176.798997  -177.149058  -177.149058  -177.822663
Ehar        31.978599      36.142743    36.329001    36.329001    36.410612
Eexc        0.0            -7.194596    -7.765807    -7.765807    -8.267519
Etot        -67.038565     -73.925425   -74.470692   -74.470692   -74.945193

OCCUPIED ORBITALS (energy / occupation)
shell   Hartree only    LDA x-only      LDA             LDA (Libxc)     PBE
1s      -16.9125 / 2.0  -18.6908 / 2.0  -18.7582 / 2.0  -18.7582 / 2.0  -18.8986 / 2.0
2s      -0.5229 / 2.0   -0.8206 / 2.0   -0.8712 / 2.0   -0.8712 / 2.0   -0.8788 / 2.0
2p      -0.0473 / 4.0   -0.2895 / 4.0   -0.3383 / 4.0   -0.3383 / 4.0   -0.3321 / 4.0
```

Reading the progression left to right:
- **Hartree only** (no XC): pure mean-field electron-electron repulsion.
  `Ehar ≈ +32 Ha`, `Etot = -67.04 Ha`. This over-binds nothing and
  over-repels everything -- there's no exchange hole keeping same-spin
  electrons apart, so every orbital sits too high (1s at `-16.9 Ha`,
  shallower than it should be) and 2s/2p aren't yet split by exchange.
- **LDA exchange only** adds Slater exchange on top of Hartree, with no
  correlation yet: `Eexc = -7.19 Ha` (exchange alone), `Etot = -73.93 Ha`.
  Exchange is the dominant piece of `Eexc` -- correlation (see next row)
  only adds another `-0.57 Ha` on top, but exchange alone already recovers
  most of the missing binding.
- **LDA** adds Perdew-Wang correlation on top of exchange: `Eexc` deepens to
  `-7.77 Ha`, `Etot = -74.47 Ha`, close to the atom's real (non-relativistic,
  all-electron) energy.
- **LDA via Libxc** uses the exact same physical functionals (Slater
  exchange + PW92 correlation) but evaluated through Libxc's compiled C
  implementation instead of `BuiltinFunctional`'s pure Julia one. Every
  reported quantity (`Etot`, all orbital energies) agrees with the `lda` row
  to the full precision shown -- confirms `BuiltinFunctional`'s `SlaterXα`/
  `PerdewWang` reproduce Libxc exactly for this atom (see also the
  dedicated unit tests in `test/physics.jl`, which check this pointwise
  against Libxc across a range of densities).
- **PBE** (GGA) makes a comparatively small further correction on this atom
  (`Etot: -74.47 -> -74.95 Ha`) -- the gradient correction sharpens
  exchange-correlation near the nucleus without changing the qualitative
  picture LDA already gives.

`density_comparison.pdf` shows the same story spatially: Hartree-only has
the most diffuse density (nothing pulling same-spin electrons apart into a
tighter exchange hole), and the three exchange-correlation cases (LDA
exchange-only, LDA, LDA-via-Libxc) sit closer together with a shorter tail.
`lda` and `lda_libxc` are visually indistinguishable in every plot, as
expected from the numbers above.

## A real bug found along the way

Building this example (the very first `lda` iteration, before any density
has been computed) surfaced a genuine bug in `PerdewWang.jl`: the
**spin-unpolarized** `eval_zk`/`eval_vrho(::PerdewWang, ρ::Real)` didn't
guard against `ρ=0`, unlike their spin-polarized counterparts
(`eval_vrho_up`/`eval_vrho_down`, which already special-case `iszero(ρup+ρdown)`).
At `ρ=0`, `rs_of_ρ(ρ) = (3/(4πρ))^(1/3) = Inf`, and `G(Inf, ...)` evaluates
to `-Inf * log(1) = NaN` -- which then poisons the assembled Hamiltonian on
the very first SCF iteration (the trial density starts at zero everywhere)
and crashes the eigensolver with "matrix contains Infs or NaNs". Any LDA run
with Perdew-Wang correlation and a zero-initialized density hit this
unconditionally, before this example ran into it.

Fixed in `src/physics/exchange correlation/PerdewWang.jl` by adding the same
`iszero(ρ) && return zero(ρ)` guard already used on the spin-polarized path
(correlation energy/potential correctly vanish as `ρ -> 0`, just like
exchange already does in `SlaterXa.jl`). Confirmed against the full test
suite (`Pkg.test()`), which passes unchanged.

## Running it

`run.jl` activates the shared `examples/` environment itself
(`Pkg.activate(joinpath(@__DIR__, ".."))`), so just run:
```julia
julia run.jl
```
