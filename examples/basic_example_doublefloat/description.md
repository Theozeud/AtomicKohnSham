# Description

`examples/basic_example` redone with all FEM/SCF arithmetic carried out in
`Double64` (~32 significant digits, via `DoubleFloats.jl`) instead of
`Float64`: sodium (`Z=11, N=11`), `Rmax=500`, `Nmesh=20`, `ordermax=10`,
`lh=2`, `nh=10`, `OptimizedAufbau(max_degen=2, tol=0.01)`,
`ODA(tinit=0.6, scftol=1e-20)`, `maxiter=100`, same output pipeline
(`log.txt`, `results.txt`, `density.pdf`, `orbitals.pdf`, `potentials.pdf`,
`convergence.pdf`, `energy_breakdown.pdf`).

**One parameter differs from `basic_example` on purpose: no exchange
functional** (`ex = ec = NoFunctional(1)` instead of Slater exchange), for a
reason worth spelling out, because it wasn't the original plan.

## Why: Gauss-Legendre quadrature is Float64-only

The first version of this example kept Slater exchange (via
`BuiltinFunctional(:lda_x)`, since `Libxc.jl`'s `Functional` is a C-library
wrapper that only accepts `Float64` — that fix is still in
`BuiltinFunctional.jl` and still correct). But the total energy came out
capped at Float64 precision regardless: tightening `scftol` from `1e-9` to
`1e-20` barely moved the last few digits of `Etot`.

The cause: `Ekin`/`Ecou`/`Ehar` are computed by **exact polynomial
integration** and are genuinely `Double64`-precise. `Eexc`, though, is
computed via `compute_exc_energy`, which integrates using quadrature nodes
from `GaussLegendre` (`src/fem/integration methods.jl`). Those nodes come
from `FastGaussQuadrature.gausslegendre(npoints)`, which has exactly one
method and always returns `Float64` — there's no generic-precision variant.
Worse, the failure is easy to miss: intermediate quantities like `Pgenx`
*do* get computed at full `Double64` precision inside `GaussLegendre`'s
constructor (confirmed by evaluating the same expression in isolation), but
the struct's type parameter is pinned to `eltype(x)` — and `x` is the
`Float64` array straight out of `gausslegendre`. The `Double64`-precise
`Pgenx` gets silently truncated back to `Float64` the moment it's stored
into the `GaussLegendre{Float64}` struct. So any model with an
exchange-correlation functional is presently ceiling-capped at Float64
precision no matter what type the rest of the pipeline uses. Fixing that for
real means a generic-precision Gauss-Legendre node/weight generator — doable
without a new dependency (`Legendre(n,T)` in `fem/legendre polynomial.jl` is
already generic, so nodes could be its roots via Newton's method), but real
numerical work, out of scope here.

So: `hartree=1` is kept (the Hartree term is exact integration too — the `F`
tensor, not quadrature — so it has no such problem), but there's no
exchange-correlation functional. `GaussLegendre` is still constructed in
`run.jl`, matching `basic_example`'s structure, but with no XC functional
it's never actually touched by the computation.

## Cross-validation

Running the *same* Hartree-only, no-XC configuration in both precisions:

|      | Float64            | Double64                          |
|------|---------------------|------------------------------------|
| Etot | -148.38279166622362 | -148.382791666223127793970294844... |

These agree to essentially the full precision `Float64` can express (~16
significant digits, right at `eps(Float64)`) — not partway through and then
diverging, like the Slater-exchange version did. That's the direct
confirmation: with the quadrature-dependent term removed, `Double64`
actually delivers ~30 correct digits, `Float64` reproduces as much of that
as its own precision allows, and the earlier, XC-including version's
mismatch really was the `GaussLegendre` cap and nothing else.

One remaining wrinkle, for anyone reading `results.txt` closely: the
"density norm error" sanity check (`∫4πr²ρdr - N`) still shows `~2.7e-9`,
not machine-`Double64`-small like the other four sanity checks. That check
independently re-integrates the density via `GaussLegendre`'s (`Float64`)
quadrature nodes purely for verification — it's a limitation of that
diagnostic itself, not of the actual SCF computation, which the other four
checks (energy balance, electron count, orthonormality, and a stopping
criterion that reaches `1.6e-30`) confirm is genuinely full-precision.

`run.jl` activates the shared `examples/` environment itself
(`Pkg.activate(joinpath(@__DIR__, ".."))`), so just run:
```julia
julia run.jl
```
