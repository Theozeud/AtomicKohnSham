```@meta
CurrentModule = AtomicKohnSham
```

# Solving for the ground state

The nonlinear eigenvalue problem is solved by a **Self-Consistent Field
(SCF)** iteration: build an effective Hamiltonian from the current density,
diagonalize it to get new orbitals, form the new density, and repeat until
convergence.

## SCF algorithm

The **Optimal Damping Algorithm (ODA)** is currently the supported SCF
scheme. It mixes successive densities along an energy-minimizing line search
rather than replacing them outright, which stabilizes convergence:

```@docs
ODA
```

## Aufbau (occupation) scheme

`ODA` needs a strategy for turning orbital energies into occupation numbers
at each iteration — the `aufbau` keyword:

```@docs
OptimizedAufbau
```

Two alternatives exist for special cases: `FrozenAufbau(n::Dict)` fixes the
occupations to user-prescribed values (ignoring orbital energies entirely),
and `SmearedAufbau(Temp)` is a thermal-smearing scheme.

## Running the SCF loop

```@docs
groundstate
```

For finer control (inspecting the solver between iterations, custom
callbacks — see [Reports, logs & plots](@ref)), the lower-level
`KSESolver`/`solve!` pair that `groundstate` wraps is also part of the
public API.

## Example

```@example manual
using AtomicKohnSham
using Libxc

model = KSEModel(Z = 11, N = 11, ex = Functional(:lda_x; n_spin = 1), ec = NoFunctional(1))
mesh = expmesh(0, 500, 60; s = 1.2)
basis = P1IntLegendreBasis(mesh; ordermax = 10)
discretization = KSEDiscretization(basis, model; lh = 2, nh = 10,
    fem_integration_method = GaussLegendre(basis, 2000))

aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-1)
alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-9)

sol = groundstate(model, discretization, alg; maxiter = 100)
println(sol)
```

Continue to [Analyzing a solution](@ref) to inspect `sol` in more detail.
