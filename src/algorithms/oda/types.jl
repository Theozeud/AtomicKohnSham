# ===================================================================
#                   OPTIMAL DAMPLING ALGORITHM
# ===================================================================
struct ODA{T <: Real, AufbauType <: Aufbau} <: SCFAlgorithm
    scftol::T                   # SCF tolerance
    tinit::T                    # Initial relaxation parameter
    frozen_t::Bool              # True if t is not optimized
    aufbau::AufbauType          # Aufbau Method used
    maxiter_ls::Int             # Maximum iteration for linesarch
    abstol_ls::T                # Absolute tolerance for line_search
    reltol_ls::T                # Relative tolerance for linesearch

    function ODA(;scftol::Real, tinit::Real, frozen_t::Bool = false,
                aufbau::Aufbau = OptimizedAufbau(), maxiter_ls::Int = 100,
                abstol_ls = 10^4*eps(typeof(tinit)), reltol_ls = abstol_ls)

        new{typeof(tinit),
            typeof(aufbau)}(scftol, tinit, frozen_t, aufbau, maxiter_ls,
                            abstol_ls, reltol_ls)
    end
end


# ===================================================================
#                               CACHE
# ===================================================================
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


function create_cache_alg(alg::ODA, discretization::KSEDiscretization, model::KSEModel)
    elT = eltype(discretization)
    t = elT(alg.tinit)
    D = zero_density(discretization)
    Dprev = zero_density(discretization)
    Dbuf = zero_density(discretization)
    U = zero_orbitals_coeffs(discretization)
    ϵ = zero_orbitals_energies(discretization)
    n = zero_occupation_numbers(discretization)
    energies_prev = Energies(elT)
    aufbaucache = OptimizedAufbauCache(discretization, model, alg.aufbau)

    F = if has_exchcorr(model)
        t -> begin
            @. Dbuf = t*D + (1-t)*Dprev
            compute_exc_energy(discretization, model, Dbuf)
        end
    else
        zero(elT)
    end

    ODACache{ elT, typeof(D), typeof(U), typeof(ϵ), typeof(n), typeof(aufbaucache),
        typeof(F)}(t, D, Dprev, Dbuf, U, ϵ, n, energies_prev, aufbaucache, F)
end

# ===================================================================
#                          ODA SOLUTION
# ===================================================================
struct ODASolution{ densityType <: AbstractArray,
                    orbitalsType <: AbstractArray,
                    orbitalsenergyType <: AbstractArray,
                    occupationType,
                    nType <: AbstractArray} <: SCFSolution
    density_coeffs::densityType                 # Density Matrix at final time
    orbitals::orbitalsType                      # Coefficient of orbitals at final time
    orbitals_energy::orbitalsenergyType         # Orbitals energy at final time
    occupation_number::occupationType           # Occupation number at final time
    n::nType
end

function makesolution(cache::ODACache, ::ODA, solver::KSESolver)
    @unpack ϵ, n = cache
    index = findall(x->x ≠ 0, n)
    index_sort = sortperm(ϵ[index])
    new_index = index[index_sort]
    occupation = if solver.discretization.nspin == 1
        [(string(i[2] + i[1] - 1, L_QUANTUM_LABELS[i[1]]), ϵ[i], n[i])
        for i in new_index]
    else
        [(string(i[2] + i[1] - 1, L_QUANTUM_LABELS[i[1]], SPIN_LABELS[i[3]]), ϵ[i], n[i])
        for i in new_index]
    end
    ODASolution{typeof(cache.D), typeof(cache.U), typeof(cache.ϵ),
    typeof(occupation), typeof(cache.n)}(
        cache.D, cache.U, cache.ϵ, occupation, cache.n)
end


# ===================================================================
#                               ODA LOG
# ===================================================================
#=
struct ODALog
    orbitals::Any
    orbitals_energy::Any
    density::Any
    occupation_number::Any
end

function create_logbook(::ODA)
    orbitals = []
    orbitals_energy = []
    density = []
    occupation_number = []
    ODALog(orbitals, orbitals_energy, density, occupation_number)
end
=#
