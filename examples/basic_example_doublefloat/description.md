# Description

`examples/basic_example` redone with all FEM/SCF arithmetic carried out in
`Double64` (~32 significant digits, via `DoubleFloats.jl`) instead of
`Float64`, with **exactly the same** physical, discretization, and algorithm
parameters as the original: sodium (`Z=11, N=11`), Slater exchange only,
`Rmax=500`, `Nmesh=20`, `ordermax=10`, `lh=2`, `nh=10`,
`GaussLegendre(basis, 2000)`, `OptimizedAufbau(max_degen=2, tol=0.1)`,
`ODA(tinit=0.6, scftol=1e-9)`, `maxiter=100`. Same output pipeline too:
`log.txt`, `results.txt`, `density.pdf`, `orbitals.pdf`, `potentials.pdf`,
`convergence.pdf`, `energy_breakdown.pdf`.

## What went into getting here

This took a few iterations, worth recording since the intermediate failures
are as informative as the final result.

**1. Libxc can't do extended precision.** `Functional(:lda_x)` from
`Libxc.jl` calls into a compiled C library that only accepts `Float64`
arrays. This example uses `BuiltinFunctional(:lda_x)` instead —
`AtomicKohnSham`'s own pure-Julia Slater exchange. Its docstring already
claimed "higher precision than Float64"; that turned out not to be true yet
(see point 3), fixed along the way.

**2. `Z`, `N` are plain `Int`.**
They used to have to be: `resolve_twofold_degeneracy!` (used here, since
3s/3p end up nearly degenerate) had a hardcoded `Nleft::Int` argument. It now
takes `Nleft::Real`, so `Double64` (or non-integer) `Z`, `N` work as well.
Keeping them as `Int` does not affect precision: the mesh, basis and all SCF
arithmetic are controlled separately via `expmesh(...; T=Double64)` and
`P1IntLegendreBasis(mesh, Double64; ...)`.

**3. The real blocker: `GaussLegendre` was silently Float64-only.**
`Ekin`/`Ecou`/`Ehar` are exact polynomial integration and were always
genuinely `Double64`-precise. `Eexc`, though, is computed via
`compute_exc_energy`, which integrates using `GaussLegendre`'s quadrature
nodes — and those came from `FastGaussQuadrature.gausslegendre(npoints)`,
which has exactly one method and always returns `Float64`, with no
generic-precision variant. The failure was subtle: intermediate quantities
like `Pgenx` were genuinely computed at `Double64` precision inside the
constructor, then silently truncated back to `Float64` the moment they were
stored into the struct, because the struct's type parameter was inferred
from `eltype(x)` — and `x` could only ever be `Float64`. First attempt at
this example (dropping the exchange functional entirely, keeping only the
exact-integration Hartree term) sidestepped the bug rather than fixing it.

**The actual fix** (in `src/fem/integration methods.jl`): `GaussLegendre`'s
constructor now dispatches on `T = eltype(basis)` — `Float64` still goes
through `FastGaussQuadrature.gausslegendre` (fast: ~O(n), unchanged for
every existing user), anything else goes through `QuadGK.gauss(T, npoints)`
(generic-precision Golub-Welsch, `O(n²)`, added as a new dependency). This
matters for speed: benchmarked at `npoints=2000`, `FastGaussQuadrature`
takes ~0.02ms, `QuadGK.gauss(Float64, 2000)` takes ~408ms — about 20,000x
slower at equal `n`, and `QuadGK.gauss(Double64, 2000)` takes ~1s. That's
why the dispatch matters rather than just always using `QuadGK`: `Float64`
users pay nothing, and the `O(n²)`/~1s cost only applies where it's
unavoidable anyway (`GaussLegendre` is built once per discretization, not
once per SCF iteration, so a one-time ~1s setup cost is negligible next to
a calculation that runs for tens of seconds to minutes).

## Cross-validation

|      | `basic_example` (Float64) | `basic_example_doublefloat` |
|------|---------------------------|------------------------------|
| Etot | -160.6282275965646        | -160.628227596564585699661445675598771 |

These now agree to essentially the full precision `Float64` can express
(~16 significant digits, right at `eps(Float64)`) — not partway through and
diverging early like the pre-fix version did. That's the confirmation the
fix is real: two independent-precision code paths computing the same
physics landing on the same answer, each accurate to its own precision
limit.

One remaining wrinkle: the "density norm error" sanity check in
`results.txt` (`∫4πr²ρdr - N`) still shows `~4e-9`, not machine-`Double64`
small. That check independently re-integrates the density using
`GaussLegendre`'s quadrature nodes purely for verification — with the fix
in place those nodes *are* genuinely `Double64`-precise now, so `~4e-9` here
reflects the quadrature rule's truncation error at `npoints=2000` for this
particular (non-polynomial) integrand, not a precision bug. The other four
sanity checks (energy balance ~3e-30, electron count exact, orthonormality
~3e-25, and a stopping criterion that reaches the requested `scftol`)
confirm the SCF computation itself is fully precise.

`run.jl` activates the shared `examples/` environment itself
(`Pkg.activate(joinpath(@__DIR__, ".."))`), so just run:
```julia
julia run.jl
```
