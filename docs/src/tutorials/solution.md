```@meta
CurrentModule = AtomicKohnSham
```

# Analyzing a solution

`groundstate` returns a [`KSESolution`](@ref), a self-contained snapshot of
the converged (or stopped) calculation: energies, orbital coefficients,
occupations, and enough context (model, discretization, basis) to evaluate
anything else on demand.

```@docs
KSESolution
```

Printing a solution (`println(sol)`, or just `sol` at the REPL) shows its
name, convergence status, and occupied orbitals with their energies and
occupation numbers — see the example at the end of
[Solving for the ground state](@ref).

## Setup

This continues the running example (each manual page is self-contained, so
we solve it again here):

```@example manual
using AtomicKohnSham
using Libxc

model = KSEModel(Z = 11, N = 11, ex = Functional(:lda_x; n_spin = 1), ec = NoFunctional(1))
mesh = expmesh(0, 500, 60; s = 1.2)
basis = P1IntLegendreBasis(mesh; ordermax = 10)
discretization = KSEDiscretization(basis, model; lh = 2, nh = 10,
    fem_integration_method = GaussLegendre(basis, 2000))
alg = ODA(tinit = 0.6, aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-1), scftol = 1e-9)

sol = groundstate(model, discretization, alg; maxiter = 100)
nothing # hide
```

## Energies

`sol.energies` is an [`Energies`](@ref) holding the kinetic, nuclear
attraction (`Ecou`), Hartree, exchange–correlation, and total energy:

```@docs
Energies
```

```@example manual
(Etot = sol.energies.Etot, Ekin = sol.energies.Ekin, Ecou = sol.energies.Ecou,
 Ehar = sol.energies.Ehar, Eexc = sol.energies.Eexc)
```

## Evaluating orbitals and the density

Radial quantities are evaluated pointwise on a vector of radii `X` (which
should not include `r = 0`, a removable singularity for the quantities
below):

```@docs
eval_orbital
eval_density
```

```@example manual
X = [0.1, 0.5, 1.0, 2.0, 5.0]
eval_orbital(sol, "3s", X)
```

```@example manual
eval_density(sol, X)
```

## Evaluating potentials

The individual contributions to the effective potential, and their sum, can
also be evaluated pointwise — useful for plotting the potential well an
orbital sits in (see [Reports, logs & plots](@ref)):

```@docs
eval_hartree
eval_nuclear
eval_kinetic_potential
eval_vxc
eval_effective_potential
```

```@example manual
eval_effective_potential(sol, 0, X)  # l = 0 channel
```

Continue to [Reports, logs & plots](@ref) to export `sol` to text and plots.
