"""
    SanityChecks{T}

Physical consistency diagnostics computed from a converged (or stopped)
Kohn–Sham solution, independent of any comparison to a previous run's output.

Each field checks an identity or conservation law that must hold regardless of
the specific system/functional/discretization. A nonzero value here indicates
an actual computational problem (in FEM assembly, the eigensolver, the aufbau
scheme, or energy bookkeeping) even when the SCF loop reports `SUCCESS` and no
exception is thrown along the way — the kind of bug a golden-value regression
test on the final energy alone cannot distinguish from "the physics is right".

# Fields
- `energy_balance::T`: `Etot - (Ekin+Ecou+Ehar+Eexc)`. Zero by definition;
  any deviation means the energy components don't actually add up to the
  reported total.
- `electron_count_error::T`: `∑ᵢnᵢ - N`. Zero: occupation numbers must sum to
  the model's electron count.
- `density_norm_error::T`: `∫4πr²ρ(r)dr - N`, evaluated by numerical
  quadrature independently of the occupation numbers. Zero: cross-checks that
  the density matrix is consistent with the occupations, through a different
  code path than `electron_count_error`.
- `orthonormality_error::T`: `maxₗ,σ ‖Uₗ,σᵀ M₀ Uₗ,σ - I‖∞`, where `M₀` is the
  FEM mass matrix. Zero: the Kohn–Sham orbitals within each angular
  momentum/spin channel must be orthonormal.
- `virial_ratio::T`: `Ekin / (-Etot)`. By the quantum virial theorem this goes
  to exactly `1` whenever the only potentials present are homogeneous of
  degree -1 in `r` — which includes *both* the nuclear Coulomb term and the
  Hartree term (any `hartree` coefficient), just not exchange–correlation
  (LDA-type functionals scale differently, e.g. as `ρ^(4/3)`, so no XC means
  no closed-form target). Unlike the other fields above, this is *not* a
  fixed-tolerance bookkeeping identity: in a finite radial domain `[0,Rmax]`
  it converges to `1` only as `Rmax → ∞` (the deviation is a domain-truncation
  effect — empirically it shrinks rapidly with `Rmax`, consistent with cutting
  off the exponentially decaying tail of a bound state, and is largely
  insensitive to basis/polynomial order at fixed `Rmax`), so treat it as a
  convergence diagnostic with a loose tolerance, not an exact check.

Any field is `NaN` if the underlying quantity (`D`, `U`, or `n`) is
unavailable from the algorithm cache used to solve the problem.
"""
struct SanityChecks{T <: Real}
    energy_balance::T
    electron_count_error::T
    density_norm_error::T
    orthonormality_error::T
    virial_ratio::T
end

function _density_norm_error(discretization::KSEDiscretization{T}, D) where {T}
    isnothing(D) && return T(NaN)
    @unpack fem_integration_method, nspin, N = discretization
    @unpack y, wy = fem_integration_method
    ρ = similar(y)
    if nspin == 1
        eval_density!(ρ, discretization, D, y)
    else
        @views DUP = D[:, :, 1]
        @views DDOWN = D[:, :, 2]
        ρup = similar(y)
        ρdown = similar(y)
        eval_density!(ρup, discretization, DUP, y)
        eval_density!(ρdown, discretization, DDOWN, y)
        @. ρ = ρup + ρdown
    end
    @. ρ *= y^2
    4π * dot(wy, ρ) - N
end

function _orthonormality_error(discretization::KSEDiscretization{T}, U) where {T}
    isnothing(U) && return T(NaN)
    @unpack lₕ, nₕ, nspin = discretization
    M₀ = discretization.femops.M₀
    err = zero(T)
    for σ in 1:nspin
        for l in 1:(lₕ + 1)
            Ulσ = nspin == 1 ? (@view U[:, :, l]) : (@view U[:, :, l, σ])
            gram = Ulσ' * M₀ * Ulσ
            err = max(err, maximum(abs, gram - I))
        end
    end
    err
end

function _sanity_checks(discretization::KSEDiscretization{T}, model::KSEModel,
        D, U, n, energies::Energies) where {T}
    energy_balance = energies.Etot -
                      (energies.Ekin + energies.Ecou + energies.Ehar + energies.Eexc)
    electron_count_error = isnothing(n) ? T(NaN) : sum(n) - model.N
    density_norm_error = _density_norm_error(discretization, D)
    orthonormality_error = _orthonormality_error(discretization, U)
    virial_ratio = iszero(energies.Etot) ? T(NaN) : energies.Ekin / (-energies.Etot)
    SanityChecks{T}(energy_balance, electron_count_error, density_norm_error,
        orthonormality_error, virial_ratio)
end
