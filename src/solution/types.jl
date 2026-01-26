"""
    KSESolution(solver::KSESolver, name::String) -> KSESolution

Constructs a solution object from a fully converged (or stopped) Kohn–Sham solver.
This structure stores the result of a self-consistent field (SCF) computation,
including metadata such as convergence status, final energies, iteration logs,
and algorithm-specific data.

# Arguments
- `solver::KSESolver`: The SCF solver after a call to `solve!`.
- `name::String`: Name or label for the solution (can be empty if unused).

# Fields
- `success::String`: Convergence status — either `"SUCCESS"` if convergence occurred before reaching `maxiter`,
                                            or   `"MAXITERS"` if the iteration limit was reached.
- `solveropts::SolverOptions`: Options that were used in the solver (e.g., tolerance, integration method).
- `niter::Int`: Number of SCF iterations performed.
- `stopping_criteria::T`: Final convergence metric value (e.g., norm of residual density).
- `energies::Dict{Symbol, T}`: Energies computed (total, kinetic, Hartree, etc.).
- `datas::SCFSolution`: Struct with detailed solution data, whose fields depend on the SCF algorithm used.
- `log::LogBook`: Object that records iteration history and diagnostics.
- `name::String`: Name of the solution (can be used for identification, labeling results, etc.).
"""
struct KSESolution{ T <: Real, TD, TU, TE, TN, TO, logbookType <: LogBook, C<:KSEContext}
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

    logbook::logbookType                # LogBook of all recorded quantities

    name::String                        # Name of the solution

    context::C                          # Minimal context on the problem solved

    function KSESolution(solver::KSESolver, name::String)

        @unpack model, alg, logbook, discretization = solver
        @unpack lh, nh, fem_integration_method

        # Status of the solver
        success = solver.niter == solver.maxiter ? "MAXITERS" : "SUCCESS"

        # Datas specific to the algorithm
        D = getfield_or_nothing(solver.algcache, :D)
        U = getfield_or_nothing(solver.algcache, :U)
        ϵ = getfield_or_nothing(solver.algcache, :ϵ)
        n = getfield_or_nothing(solver.algcache, :n)
        occupied = occupied_orbitals_summary(discretization, n, ϵ)

        # Context
        context = KSEContext(model, alg, lh, nh, basis, fem_integration_method)

        new{eltype(discretization),
            typeof(D),
            tyepof(U),
            typeof(ϵ),
            typeof(n),
            typeof(occupied),
            typeof(logbook),
            typeof(context)}(success,
                             solver.niter,
                             solver.stopping_criteria,
                             solver.energies,
                             D, U, ϵ, n, occupied,
                             solver.logbook,
                             name,
                             context)
    end
end

"""
Return the Kohn–Sham orbital coefficients stored in the solution.
"""
orbitals(sol) = sol.U

"""
Return the Kohn–Sham orbital energies.
"""
orbital_energies(sol) = sol.ϵ

"""
Return the orbital occupation numbers.
"""
occupations(sol) = sol.n

"""
Return the coefficients of the one-particle reduced density.
"""
density_matrix(sol) = sol.D


@generated function getfield_or_nothing(x, ::Val{s}) where {s}
    if s ∈ fieldnames(x)
        :(getfield(x, s))
    else
        :(nothing)
    end
end


function occupied_orbitals_summary(discretization::KSEDiscretization, n::AbstractArray{T},
                                  ϵ::AbstractArray{T}) where T <: Real
    index = findall(x->x ≠ 0, n)
    index_sort = sortperm(ϵ[index])
    new_index = index[index_sort]
    [(shell_string(convert_index_nl(discretization, i)), ϵ[i], n[i]) for i in new_index]
end

function occupied_orbitals_summary(::KSEDiscretization, ::Nothing, ::Nothing)
    nothing
end
