```@meta
CurrentModule = AtomicKohnSham
```

# Discretization

Assuming spherical symmetry, each Kohn–Sham orbital factors into a radial
part `u(r)` and a spherical harmonic. `AtomicKohnSham.jl` discretizes the
radial part `u(r)` with a **finite element method (FEM)**: a mesh cuts the
radial domain `[0, Rmax]` into cells, and a polynomial basis is built on
top of it.

## Mesh

Several radial mesh types are available; `expmesh` (exponentially graded,
denser near the origin where orbitals vary fastest) is recommended for
atomic problems:

```@docs
Mesh
expmesh
```

`linmesh`, `geometricmesh`, `polynomialmesh` and `explinmesh` follow the same
`(a, b, n; kwargs...) -> Mesh` convention for, respectively, a uniform,
geometric, polynomially-graded, or exponential-then-linear mesh.

## FEM basis

Currently, the only implemented basis is a hierarchical **integrated
Legendre polynomial** basis (`P1` vertex functions plus higher-order
"bubble" functions per cell, up to `ordermax`):

```@docs
P1IntLegendreBasis
```

## Integration method

FEM matrix assembly and energy integrals need numerical quadrature:

```@docs
GaussLegendre
```

## Putting it together

```@docs
KSEDiscretization
```

`lh` (the angular momentum cutoff) should cover every occupied shell's `l`
— e.g. `lh = 1` for an atom occupying up to a `p` shell. `nh` bounds how
many radial solutions per channel are kept.

## Example

This continues the running example from [The physical model](@ref) (each
manual page is self-contained, so we redefine `model` here too):

```@example manual
using AtomicKohnSham
using Libxc

model = KSEModel(Z = 11, N = 11, ex = Functional(:lda_x; n_spin = 1), ec = NoFunctional(1))

mesh = expmesh(0, 500, 60; s = 1.2)
basis = P1IntLegendreBasis(mesh; ordermax = 10)
discretization = KSEDiscretization(basis, model; lh = 2, nh = 10,
    fem_integration_method = GaussLegendre(basis, 2000))
discretization.Nₕ
```

Continue to [Solving for the ground state](@ref) to run the SCF loop.
