```@meta
CurrentModule = AtomicKohnSham
```

# The physical model

A [`KSEModel`](@ref) collects everything that defines *what* problem is being
solved: the nucleus, the number of electrons, and which exchange–correlation
approximation stands in for the true electron–electron interaction beyond
the (explicit) Hartree term.

```@docs
KSEModel
```

## Exchange–correlation functionals

`ex` and `ec` accept any of:

- **[Libxc.jl](https://github.com/unkcpz/Libxc.jl) functionals**, e.g.
  `Functional(:lda_x; n_spin = 1)` for Slater exchange, or
  `Functional(:lda_c_pw; n_spin = 2)` for Perdew–Wang correlation in a
  spin-polarized (LSDA) calculation. This gives access to the full LDA/LSDA
  catalog implemented by Libxc.
- **`NoFunctional(n_spin)`**, to disable exchange or correlation entirely
  (e.g. for a reduced Hartree–Fock-like model).
- **Built-in functionals** (`BuiltinFunctional`), implemented directly in
  this package to be run at arbitrary precision (`Double64`, etc.), which
  Libxc's `Float64`-only C implementation cannot do. Currently only Slater
  exchange (`:lda_x`) is available this way.

`ex`/`ec` must agree on the number of spin channels; `model.nspin` is
inferred automatically as `max(ex.n_spin, ec.n_spin)`.

## Convenience constructors

For the two most common setups:

```@docs
RHF
Slater
```

## Example

The rest of this manual builds up a single running example — the ground
state of a spin-unpolarized sodium atom with the Slater exchange functional
— one piece per page.

```@example manual
using AtomicKohnSham
using Libxc

model = KSEModel(Z = 11, N = 11, ex = Functional(:lda_x; n_spin = 1), ec = NoFunctional(1))
model.nspin
```

Continue to [Discretization](@ref) to turn this model into a solvable
finite-dimensional problem.
