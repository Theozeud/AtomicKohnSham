"""
Return `getfield(x, s)`, or `nothing` if `x` has no field `s`.

Used to pull algorithm-specific quantities (e.g. `:D`, `:U`, `:ϵ`, `:n`, `:t`)
out of an `SCFCache` without each algorithm having to implement a common
accessor interface.
"""
getfield_or_nothing(x, s::Symbol) = getfield_or_nothing(x, Val(s))
@generated function getfield_or_nothing(x, ::Val{s}) where {s}
    if s ∈ fieldnames(x)
        :(getfield(x, s))
    else
        :(nothing)
    end
end

"""
Summarize the occupied orbitals of a discretization as
`(shell_label, orbital_energy, occupation_number)` tuples, sorted by orbital
energy. Returns `nothing` if occupations/energies are unavailable (e.g. the
algorithm cache has no `:n`/`:ϵ` fields).
"""
function _occupied_orbitals_summary(discretization::KSEDiscretization, n::AbstractArray{T},
                                  ϵ::AbstractArray{T}) where T <: Real
    index = findall(x->x ≠ 0, vec(n))
    index_sort = sortperm(ϵ[index])
    new_index = index[index_sort]
    [(shell_string(convert_index_nl(discretization, i)), ϵ[i], n[i]) for i in new_index]
end

function _occupied_orbitals_summary(::KSEDiscretization, ::Nothing, ::Nothing)
    nothing
end

"""
    KSESolution(solver::KSESolver, name::String) -> KSESolution

Construct a solution object from a Kohn–Sham solver after a call to `solve!`.

This structure stores the final state of a self-consistent field (SCF)
computation, including convergence information, final energies, density and
orbital coefficients, orbital energies, occupation numbers, Hartree potential
coefficients, iteration logs, and the minimal context needed to post-process the
solution.

# Arguments
- `solver::KSESolver`: Solver containing the final state of the SCF computation.
- `name::String`: Name or label associated with the solution.

# Fields
- `success::String`: Final solver status. It is `"SUCCESS"` if the solver stopped
  before reaching `maxiter`, and `"MAXITERS"` if the maximum number of iterations
  was reached.
- `niter::Int`: Number of SCF iterations performed.
- `stopping_criteria::T`: Final value of the stopping criterion.
- `energies::Energies{T}`: Final energy contributions stored in an `Energies`
  object.
- `D`: Coefficients of the one-body reduced density matrix.
- `U`: Coefficients of the Kohn–Sham orbitals in the FEM basis.
- `ϵ`: Kohn–Sham orbital energies.
- `n`: Orbital occupation numbers.
- `occupied`: Summary of occupied orbitals, sorted by orbital energy. Each entry
  contains the shell label, orbital energy, and occupation number.
- `W`: Coefficients of the Hartree potential.
- `sanity::SanityChecks`: Physical consistency diagnostics (energy bookkeeping,
  electron count, density normalization, orbital orthonormality, virial
  ratio) computed independently of any comparison to a previous run's output.
- `logbook::LogBook`: Logbook containing recorded quantities during the SCF
  iterations.
- `name::String`: Name of the solution, used for identification or labeling.
- `context::KSEContext`: Minimal context of the solved problem, containing the
  model, algorithm, discretization sizes, FEM basis, and integration method.
"""
struct KSESolution{T <: Real, TD, TU, TE, TN, TO, TW, logbookType <: LogBook, C<:KSEContext}
    success::String                     # Print the final state of the solver
                                        # Can be : SUCCESS or MAXITERS

    niter::Int                          # Number of iterations
    stopping_criteria::T                # Final stopping criteria

    energies::Energies{T}               # Energies

    D::TD                               # Coefficients of the one-body reduced density
    U::TU                               # Coefficients of the orbitals
    ϵ::TE                               # Energies of the orbitals
    n::TN                               # Occupations numbers
    occupied::TO                        # Occupied orbitals with energies

    W::TW                               # Coefficients of the Hartree potential

    sanity::SanityChecks{T}             # Physical consistency diagnostics

    logbook::logbookType                # LogBook of all recorded quantities

    name::String                        # Name of the solution

    context::C                          # Minimal context on the problem solved

    function KSESolution(solver::KSESolver, name::String)

        @unpack model, alg, logbook, discretization = solver
        @unpack lₕ, nₕ, fem_integration_method, basis = discretization

        # Status of the solver
        success = solver.niter == solver.maxiter ? "MAXITERS" : "SUCCESS"

        # Datas specific to the algorithm
        D = getfield_or_nothing(solver.algcache, :D)
        U = getfield_or_nothing(solver.algcache, :U)
        ϵ = getfield_or_nothing(solver.algcache, :ϵ)
        n = getfield_or_nothing(solver.algcache, :n)
        occupied = _occupied_orbitals_summary(discretization, n, ϵ)

        W = discretization.cache.hartw.W

        sanity = _sanity_checks(discretization, model, D, U, n, solver.energies)

        # Context
        context = KSEContext(model, alg, lₕ, nₕ, basis, fem_integration_method)

        new{eltype(discretization),
            typeof(D),
            typeof(U),
            typeof(ϵ),
            typeof(n),
            typeof(occupied),
            typeof(W),
            typeof(logbook),
            typeof(context)}(success,
                             solver.niter,
                             solver.stopping_criteria,
                             solver.energies,
                             D, U, ϵ, n, occupied,
                             W,
                             sanity,
                             solver.logbook,
                             name,
                             context)
    end
end
