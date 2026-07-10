```@meta
CurrentModule = AtomicKohnSham
```

# AtomicKohnSham.jl

**AtomicKohnSham.jl** computes the ground state of isolated atoms and ions
within the framework of **Extended Kohn–Sham models**. It was originally
built to investigate the existence of negative ions, and doubles as a
testbed for new density functionals and for exploring numerical precision
questions in atomic structure theory.

## How it works

Extended Kohn–Sham models turn the many-body Schrödinger equation into a
nonlinear eigenvalue problem: a set of one-particle orbitals, solved
self-consistently against an effective potential that itself depends on the
orbitals. Assuming spherical symmetry (true for an isolated atom or ion),
that 3D problem reduces to a family of 1D radial equations, one per angular
momentum channel.

This package solves those radial equations with a **high-precision finite
element method** (piecewise polynomials of degree up to 20 on graded
meshes), which reaches high accuracy with comparatively few mesh points —
important for resolving open numerical questions, such as whether the
energy of the outermost orbital of certain ions is exactly zero or merely
very close to it.

## Installation

```julia
] add AtomicKohnSham
```

## A first calculation

Every calculation follows the same four steps: describe the physical
**model**, choose a **discretization**, pick an SCF **algorithm**, then
**solve**. Here is the ground state of hydrogen with the Slater exchange
functional (no correlation):

```@example quickstart
using AtomicKohnSham
using Libxc

# 1. Model: nuclear charge, electron count, exchange/correlation functionals
model = KSEModel(Z = 1, N = 1, ex = Functional(:lda_x; n_spin = 1), ec = NoFunctional(1))

# 2. Discretization: radial mesh + FEM basis
mesh = expmesh(0, 30, 20; s = 1.5)
basis = P1IntLegendreBasis(mesh; ordermax = 10)
discretization = KSEDiscretization(basis, model; lh = 0, nh = 3)

# 3. Algorithm: Optimal Damping Algorithm (ODA)
alg = ODA(tinit = 0.6, scftol = 1e-10)

# 4. Solve
sol = groundstate(model, discretization, alg; maxiter = 100)
println(sol)
```

The rest of this manual walks through each of these four steps in more
detail, then covers how to inspect a solution (energies, orbitals, density,
potentials) and how to export it (text reports, iteration logs, and plots):

1. [The physical model](@ref)
2. [Discretization](@ref)
3. [Solving for the ground state](@ref)
4. [Analyzing a solution](@ref)
5. [Reports, logs & plots](@ref)

or jump straight to the [API Reference](@ref) for the full function index.
