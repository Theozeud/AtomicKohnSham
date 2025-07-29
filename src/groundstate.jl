#--------------------------------------------------------------------
#                           GROUNDSTATE
#--------------------------------------------------------------------
"""
    groundstate(model::KSEModel,
                discretization::KSEDiscretization,
                method::SCFAlgorithm;
                name::String = "", kwargs...) -> KSESolution

Computes the ground state solution of a Kohn–Sham model using a specified discretization and SCF method.

This function creates a `KSESolver` object with the given model, discretization, and method, then solves the nonlinear eigenvalue problem via `solve!`. The result is returned as a `KSESolution` object.

# Arguments
- `model::KSEModel`: The Kohn–Sham extended model.
- `discretization::KSEDiscretization`: Discretization of the radial Kohn–Sham equations.
- `method::SCFAlgorithm`: Self-consistent field (SCF) method to solve the nonlinear problem (e.g., `ODA()`, `Quadratic()`).
- `name::String` (optional): Name for the simulation, stored in the output.
- `kwargs`: Additional keyword arguments forwarded to `KSESolver`.

# Returns
- `KSESolution`: The computed ground state, including orbitals, energies, density, etc.
"""
function groundstate(model::KSEModel,
        discretization::KSEDiscretization,
        alg::SCFAlgorithm;
        name::String = "", kwargs...)
    solver = KSESolver(model, discretization, alg; kwargs...)
    solve!(solver)
    KSESolution(solver, name)
end

"""
    groundstate(prob::AtomProblem) -> KSESolution

Convenience function to compute the ground state from an `AtomProblem`.

This function constructs the mesh, basis, and discretization from the problem specification, then calls `groundstate` with the appropriate arguments.

# Arguments
- `prob::AtomProblem`: A fully specified atomic problem including model, mesh and basis types, SCF method, cutoffs, and solver options.

# Returns
- `KSESolution`: The computed ground state solution for the atomic system.
"""
function groundstate(prob::AtomProblem)
    T = datatype(prob)
    typemesh = prob.typemesh
    typebasis = prob.typebasis
    mesh = typemesh(zero(T), prob.Rmax, prob.Nmesh; T = T, prob.optsmesh...)
    basis = typebasis(mesh, T; prob.optsbasis...)
    discretization = KSEDiscretization(prob.lh, basis, mesh, prob.model.n_spin, prob.nh)
    groundstate(prob.model, discretization, prob.alg; name = prob.name, prob.solveropts...)
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
