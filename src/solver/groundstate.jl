"""
    groundstate(model::KSEModel, discretization::KSEDiscretization, method::SCFAlgorithm;
               kwargs...) -> KSESolution

Compute the ground-state solution of a Kohn–Sham model using a given discretization and SCF
algorithm.

# Arguments
- `model::KSEModel`: Kohn–Sham extended model.
- `discretization::KSEDiscretization`: Discretization of the radial Kohn–Sham equations.
- `method::SCFAlgorithm`: SCF algorithm.

# Keyword arguments
- `maxiter::Int = 100`: Maximum number of SCF iterations.
- `callback = nothing`: Callback (or [`CallbackSet`](@ref)) invoked once per
  SCF iteration, e.g. a [`LogFileCallback`](@ref) for per-iteration diagnostics.

# Returns
- `KSESolution`: Ground-state solution including orbitals, energies, density, etc.
"""
function groundstate(model::KSEModel, discretization::KSEDiscretization, alg::SCFAlgorithm;
                    kwargs...)
    @unpack Z,N = model
    name = if  floor(Z) == Z
        charge_diff = Z - N
        if iszero(charge_diff)
            ATOMIC_NUMBER_TO_NAME[Int(Z)]
        else
            symbol_ion = charge_diff > 0 ? "+" : "-"
            ATOMIC_NUMBER_TO_NAME[Int(Z)]*string(abs(charge_diff))*symbol_ion
        end
    else
        "Z=$Z, N=$N"
    end
    solver = KSESolver(model, discretization, alg; kwargs...)
    solve!(solver)
    KSESolution(solver, name)
end


"""
    solve!(solver::KSESolver)

Run the SCF loop in-place.

At each iteration, this routine performs one SCF step (`scf_step!`), records
convergence data (`register!`), and executes the user callback (see
`callback` in [`KSESolver`](@ref), e.g. a [`LogFileCallback`](@ref) for
per-iteration diagnostics). The loop stops when the algorithm-specific
convergence test is satisfied or when `maxiter` is reached. Finally,
`postcomputations!` is called to compute any additional quantities.
"""
function solve!(solver::KSESolver)
    while solver.niter < solver.maxiter
        solver.stopping_criteria = scf_step!(solver)
        register!(solver)
        callback!(solver.callback, solver)
        solver.niter += 1
        scf_converged(solver.alg, solver) && break
    end
    postcomputations!(solver)
end


"""
    scf_step!(solver::KSESolver) -> stopping_criteria

Perform one SCF iteration and return the current stopping criterion.
"""
scf_step!(solver::KSESolver) = scf_step!(solver.alg, solver)


"""
 postcomputations!(solver::KSESolver)

Compute algorithm-specific quantities after the SCF loop.

By default this dispatches to `postcomputations!(solver.algcache, solver)` and
does nothing if the cache type does not implement postprocessing.
"""
postcomputations!(::SCFCache, ::KSESolver) = nothing
postcomputations!(solver::KSESolver) = postcomputations!(solver.algcache, solver)


"""
    register!(solver::KSESolver)

Append convergence information to the solver logbook.

The default implementation records the stopping criterion and total energy at
each iteration. Extend or replace this function to store additional data.
"""
function register!(solver::KSESolver)
    @unpack logbook = solver
    push!(logbook.stopping_criteria, solver.stopping_criteria)
    push!(logbook.Etot, solver.energies.Etot)
end
