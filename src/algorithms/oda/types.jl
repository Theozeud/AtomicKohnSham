"""
    ODA(; scftol, tinit, frozen_t=false, aufbau=OptimizedAufbau(),
         maxiter_ls=100, abstol_ls=…, reltol_ls=…)

Optimal Damping Algorithm (ODA) for self-consistent field (SCF) iterations.

`ODA` stabilizes SCF convergence by mixing successive densities using an
energy-based line search along a linear interpolation path. At each iteration,
the optimal mixing parameter is determined by minimizing a one-dimensional
energy model including kinetic, Coulomb, and Hartree contributions, with an
optional nonlinear correction (e.g. exchange–correlation).

The algorithm supports fractional occupations through a configurable Aufbau
strategy and can optionally freeze the mixing parameter.

# Keywords
- `scftol`: Convergence tolerance for the SCF iteration.
- `tinit`: Initial mixing (relaxation) parameter.
- `frozen_t`: If `true`, the mixing parameter is kept fixed to `tinit`.
- `aufbau`: Aufbau method used to fill orbital occupations.
- `maxiter_ls`: Maximum number of iterations in the line search.
- `abstol_ls`: Absolute tolerance for the line-search minimization.
- `reltol_ls`: Relative tolerance for the line-search minimization.

`ODA` is designed to be robust in the presence of degeneracies and large density
variations, and can be combined with acceleration techniques such as DIIS once
the SCF iterations are stabilized.
"""
struct ODA{T <: Real, AufbauType <: Aufbau} <: SCFAlgorithm
    scftol::T
    tinit::T
    frozen_t::Bool
    aufbau::AufbauType
    maxiter_ls::Int
    abstol_ls::T
    reltol_ls::T

    function ODA(;scftol::Real, tinit::Real, frozen_t::Bool = false,
                aufbau::Aufbau = OptimizedAufbau(), maxiter_ls::Int = 100,
                abstol_ls = 10^4*eps(typeof(tinit)), reltol_ls = abstol_ls)

        new{typeof(tinit), typeof(aufbau)}(scftol, tinit, frozen_t, aufbau, maxiter_ls,
                                           abstol_ls, reltol_ls)
    end
end


"""
    ODACache

Internal workspace for the `ODA` SCF algorithm.

Stores the current/previous densities, orbitals, occupations, and auxiliary buffers
used during ODA iterations (including line-search state and Aufbau-related data).
This cache is mutated in-place to avoid allocations inside the SCF loop.
"""

mutable struct ODACache{T <: Real,
                        DT <: AbstractArray,
                        UT <: AbstractArray,
                        ET <: AbstractArray,
                        nT <: AbstractArray,
                        AT <: AufbauCache,
                        TF} <: SCFCache

    t::T                                        # Relaxation Parameter

    D::DT                                       # Density Matrix at current time
    Dprev::DT                                   # Density Matrix at previous time
    Dbuf::DT                                    # Storage for extra Density Matrix
                                                # (usefull to avoid allocation in ODA)

    U::UT                                       # Coefficient of orbitals at current time
    ϵ::ET                                       # Orbitals energy at current time
    n::nT                                       # Occupation number at current time

    energies_prev::Energies{T}                  # Energies at previous time

    aufbaucache::AT                             # Aufbau Cache Method

    F::TF                                       # Nonlinear functional in line_search
end

function create_cache_alg(alg::ODA, discretization::KSEDiscretization{T},
                         model::KSEModel) where T
    t = T(alg.tinit)
    D = zero_density(discretization)
    Dprev = zero_density(discretization)
    Dbuf = zero_density(discretization)
    U = zero_orbitals_coeffs(discretization)
    ϵ = zero_orbitals_energies(discretization)
    n = zero_occupation_numbers(discretization)
    energies_prev = Energies(T)
    aufbaucache = create_cache_aufbau(discretization, model, alg.aufbau)

    F = if has_exchcorr(model)
        t -> begin
            @. Dbuf = t*D + (1-t)*Dprev
            compute_exc_energy(discretization, model, Dbuf)
        end
    else
        zero(T)
    end

    ODACache{T, typeof(D), typeof(U), typeof(ϵ), typeof(n), typeof(aufbaucache),
        typeof(F)}(t, D, Dprev, Dbuf, U, ϵ, n, energies_prev, aufbaucache, F)
end
