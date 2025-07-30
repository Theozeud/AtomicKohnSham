#--------------------------------------------------------------------
#                           GROUNDSTATE
#--------------------------------------------------------------------
"""
    groundstate(model::KSEModel, discretization::KSEDiscretization, method::SCFAlgorithm;
                name::String = "", kwargs...) -> KSESolution

Compute the ground-state solution of a Kohn–Sham model using a given discretization and SCF
algorithm. The result is returned as a `KSESolution`.

# Arguments
- `model::KSEModel`: Kohn–Sham extended model.
- `discretization::KSEDiscretization`: Discretization of the radial Kohn–Sham equations.
- `method::SCFAlgorithm`: SCF algorithm to solve the nonlinear problem (e.g., `ODA()`,
  `Quadratic()`).

# Keywords
- `name::String`: Optional simulation name stored in the result (default: `""`).
- `kwargs`: Additional keyword arguments forwarded to `KSESolver`.

# Returns
- `KSESolution`: Ground-state solution including orbitals, energies, density, etc.
"""
function groundstate(model::KSEModel, discretization::KSEDiscretization, alg::SCFAlgorithm;
                     name::String = "", kwargs...)
    solver = KSESolver(model, discretization, alg; kwargs...)
    solve!(solver)
    KSESolution(solver, name)
end

"""
    groundstate(prob::AtomProblem) -> KSESolution

Compute the ground-state solution of an atomic system described by an `AtomProblem`.

This convenience method builds the mesh, basis, integration method, and discretization from
the fields of the `AtomProblem`, then calls `groundstate` using the SCF algorithm and model
contained in the problem.

# Arguments
- `prob::AtomProblem`: Fully specified atomic problem, including model, mesh and basis
    types, SCF method, cutoffs, and solver options.

# Returns
- `KSESolution`: Ground-state solution for the atomic or ionic system.
"""
function groundstate(prob::AtomProblem)
    T = _datatype(prob)
    mesh = prob.typemesh(zero(T), prob.Rmax, prob.Nmesh; T = T, prob.optsmesh...)
    basis = prob.typebasis(mesh, T; prob.optsbasis...)
    model = KSEModel(;z=prob.z, N=prob.N, hartree=prob.hartree, ex=prob.ex, ec=prob.ec)
    fem_integration_method = prob.integration_method(basis, prob.optsintegration...)
    discretization = KSEDiscretization(prob.lh, basis, mesh, model.n_spin, prob.nh,
                                       fem_integration_method)
    groundstate(model, discretization, prob.alg; name = prob.name, prob.solveropts...)
end

#--------------------------------------------------------------------
#                           GROUNDSTATE STEPS
#--------------------------------------------------------------------
##
# QUESTION :
#               solver.cache et solver.alg pour dispatcher ??
#                method devrait être suffisant
##

#  SOLVE
function solve!(solver::KSESolver)
    while (solver.stopping_criteria > solver.opts.scftol || iszero(solver.niter)) &&
        solver.niter < solver.opts.maxiter
        loopheader!(solver)
        performstep!(solver.cache, solver.alg, solver)
        loopfooter!(solver)
        monitor(solver)
    end
end

# LOOPHEADER
function loopheader!(solver::KSESolver)
    # LOOPHEADER SPECIFIC FOR THE METHOD
    loopheader!(solver.cache, solver.alg, solver)
    nothing
end

# LOOPFOOTER
function loopfooter!(solver::KSESolver)
    # INCREASE THE NUMBER OF ITERATIONS DONE
    solver.niter += 1

    # COMPUTE THE NEW STOPPING CRITERIA
    stopping_criteria!(solver)

    # UPDATE THE LOG
    register!(solver)

    # LOOPFOOTER SPECIFIC FOR THE METHOD
    loopfooter!(solver.cache, solver.alg, solver)
    nothing
end

# MONITOR : DISPLAY CURRENT STATE OF SOLVER
function monitor(solver::KSESolver)
    if solver.opts.verbose > 0
        println("--------------------------")
        println("Iteration : $(solver.niter)")
    end
    if solver.opts.verbose > 1
        println("Selected Method : $(name(solver.alg))")
        println("Stopping criteria: $(solver.stopping_criteria)")
        println("Total Energy: $(solver.energies[:Etot])")
    end
    if solver.opts.verbose > 2
        monitor(solver.cache, solver.alg, solver)
    end
end

# COMPUTE THE STOPPING CRITERIA
function stopping_criteria!(solver::KSESolver)
    solver.stopping_criteria = stopping_criteria!(solver.cache, solver.alg, solver)
end

# UPDATE THE LOG : STORE INTERMEDIATE STATE
function register!(solver::KSESolver)
    @unpack logbook = solver
    @unpack stopping_criteria, energies = logbook.config

    # STORE THE STOPPING CRITERIA
    !stopping_criteria || push!(solver.logbook.stopping_criteria, solver.stopping_criteria)

    # STORE THE TOTAL ENERGY
    if energies
        for k in keys(solver.energies)
            push!(solver.logbook.energies[k], solver.energies[k])
        end
    end

    # REGISTER DATA SPECIFIC TO THE METHODS
    register!(solver.cache, solver.alg, solver)
end
