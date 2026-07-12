# Exhaustive, hand-instrumented profile of a full `groundstate` computation,
# using TimerOutputs.jl. Rather than treating `solve!`/`scf_step!` as a black
# box, this script re-implements the ODA SCF step by hand (same math, same
# function calls, same order), wrapping each conceptual phase in `@timeit` so
# the printed report shows exactly where time and memory go:
#
#   Algorithm        - bookkeeping, Aufbau occupations, ODA line search, convergence check
#   Matrix assembly  - Hartree matrix, exchange-correlation matrix, Hamiltonian sum
#   Eigenvalue problem - one dense eigensolve per (l, sigma) channel
#   Energy           - density matrix assembly + the five energy contributions
#
# This is a measurement, not a change: nothing here modifies package source.
# Run with: julia benchmark/timeroutputs_profile.jl
# Output is written to benchmark/timeroutputs_profile_output.txt

cd(@__DIR__)
using Pkg
Pkg.activate(".")

using AtomicKohnSham
using Libxc
using LinearAlgebra
using TimerOutputs

import AtomicKohnSham: KSESolver, KSESolution, register!, callback!, scf_converged,
    postcomputations!, assemble_hartree!, assemble_exc!, has_exchcorr, aufbau!,
    density!, compute_total_energy, compute_kinetic_energy, compute_coulomb_energy,
    compute_hartree_energy, compute_exc_energy, line_search_density!

"""
    timed_scf_step!(to, alg, solver)

Hand-unrolled equivalent of `AtomicKohnSham.scf_step!(alg::ODA, solver)`
(see `src/algorithms/oda/loop.jl`), with every phase wrapped in `@timeit`.
Returns the same stopping criterion as the real `scf_step!`.
"""
function timed_scf_step!(to::TimerOutput, alg::ODA, solver::KSESolver)
    discretization = solver.discretization
    model = solver.model
    energies = solver.energies
    algcache = solver.algcache
    niter = solver.niter
    D, Dprev, U, ϵ, n = algcache.D, algcache.Dprev, algcache.U, algcache.ϵ, algcache.n
    aufbaucache, energies_prev = algcache.aufbaucache, algcache.energies_prev
    ksham = discretization.ksham

    @timeit to "Algorithm" begin
        @timeit to "Bookkeeping (save previous density/energies)" begin
            @. Dprev = D
            energies_prev.Etot = energies.Etot
            energies_prev.Ekin = energies.Ekin
            energies_prev.Ecou = energies.Ecou
            energies_prev.Ehar = energies.Ehar
            energies_prev.Eexc = energies.Eexc
        end
    end

    @timeit to "Matrix assembly" begin
        @timeit to "Hartree matrix (solves A\\B)" begin
            iszero(model.hartree) ||
                assemble_hartree!(discretization, Dprev, model.N, model.hartree)
        end
        @timeit to "Exchange-correlation matrix" begin
            has_exchcorr(model) && assemble_exc!(discretization, model, Dprev)
        end
        @timeit to "Hamiltonian sum (Hfix + Vxc + Hartree)" begin
            for l in 0:discretization.lₕ
                @views vH = ksham.H[:, :, l + 1, 1]
                @. vH = ksham.Hfix[l + 1] + ksham.VxcUP + ksham.Hartree
            end
            if discretization.nspin == 2
                for l in 0:discretization.lₕ
                    @views vH = ksham.H[:, :, l + 1, 2]
                    @. vH = ksham.Hfix[l + 1] + ksham.VxcDOWN + ksham.Hartree
                end
            end
        end
    end

    @timeit to "Eigenvalue problem" begin
        for σ in 1:discretization.nspin
            for l in 0:discretization.lₕ
                @timeit to "channel (l=$l, σ=$σ)" begin
                    @views vH = ksham.H[:, :, l + 1, σ]
                    λ, V = eigen(Symmetric(discretization.femops.S * vH * discretization.femops.S))
                    @views ϵ[l + 1, :, σ] = λ[1:discretization.nₕ]
                    @views U[:, :, l + 1, σ] = discretization.femops.S * V[:, 1:discretization.nₕ]
                end
            end
        end
    end

    @timeit to "Algorithm" begin
        @timeit to "Aufbau (occupation numbers)" begin
            aufbau!(n, ϵ, U, model, discretization, aufbaucache, niter)
        end
    end

    @timeit to "Energy" begin
        if aufbaucache.postcomputations
            @timeit to "Copy trial density/energies (resolved degeneracy)" begin
                D .= aufbaucache.D1
                energies.Etot = aufbaucache.energies.Etot
                energies.Ekin = aufbaucache.energies.Ekin
                energies.Ecou = aufbaucache.energies.Ecou
                energies.Ehar = aufbaucache.energies.Ehar
                energies.Eexc = aufbaucache.energies.Eexc
            end
        else
            @timeit to "Density matrix assembly" begin
                density!(discretization, U, n, D)
            end
            @timeit to "Total energy" begin
                energies.Etot = compute_total_energy(discretization, model, D, n, ϵ)
            end
            @timeit to "Kinetic energy" begin
                energies.Ekin = compute_kinetic_energy(discretization, U, n)
            end
            @timeit to "Coulomb (nuclear attraction) energy" begin
                energies.Ecou = compute_coulomb_energy(discretization, U, n)
            end
            @timeit to "Hartree energy (solves A\\B)" begin
                energies.Ehar = compute_hartree_energy(discretization, D)
            end
            @timeit to "Exchange-correlation energy" begin
                energies.Eexc = compute_exc_energy(discretization, model, D)
            end
        end
    end

    if niter > 0
        @timeit to "Algorithm" begin
            @timeit to "Line search (ODA relaxation)" begin
                line_search_density!(alg, solver)
            end
        end
    end

    local stop
    @timeit to "Algorithm" begin
        @timeit to "Convergence check" begin
            stop_D = norm(D - Dprev)
            stop_Etot = abs(energies.Etot - energies_prev.Etot)
            stop = max(stop_D, stop_Etot)
        end
    end
    return stop
end

"""
    timed_groundstate(to, model, discretization, alg; maxiter, name)

Hand-unrolled equivalent of `AtomicKohnSham.groundstate`, timing setup, the
full SCF loop (via `timed_scf_step!`), postcomputations, and the final
`KSESolution` construction.
"""
function timed_groundstate(to::TimerOutput, model::KSEModel, discretization::KSEDiscretization,
                          alg::ODA; maxiter::Int = 100, name::String = "Benchmark")
    solver = @timeit to "Setup (init_cache! + algorithm cache)" begin
        KSESolver(model, discretization, alg; maxiter = maxiter)
    end

    @timeit to "SCF loop (total, all iterations)" begin
        while solver.niter < solver.maxiter
            solver.stopping_criteria = timed_scf_step!(to, alg, solver)
            @timeit to "Algorithm" begin
                @timeit to "Logbook + callback" begin
                    register!(solver)
                    callback!(solver.callback, solver)
                end
            end
            solver.niter += 1
            scf_converged(solver.alg, solver) && break
        end
    end

    @timeit to "Postcomputations (final Hartree potential)" begin
        postcomputations!(solver)
    end

    return @timeit to "Build KSESolution" begin
        KSESolution(solver, name)
    end
end

# ===================================================================
#                              RUN
# ===================================================================
function profile_case(Nmesh::Int; Z = 11, N = 11, Rmax = 500, ordermax = 10, lh = 2)
    model = KSEModel(Z = Z, N = N, ex = Functional(:lda_x; n_spin = 1), ec = NoFunctional(1))
    mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
    discretization = KSEDiscretization(basis, model; lh = lh, nh = 10,
        fem_integration_method = GaussLegendre(basis, 2000))
    aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-2)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-9)

    # Warm-up run (throwaway TimerOutput): forces JIT compilation of every
    # code path so the timed run below measures actual algorithmic cost, not
    # compilation latency. KSESolver re-initializes all mutable state, so
    # reusing the same discretization/model/alg for both runs is safe.
    timed_groundstate(TimerOutput(), model, discretization, alg; maxiter = 100)

    # Created only now: its internal clock must start right before the timed
    # run, otherwise setup/warm-up time leaks into "% measured" below.
    to = TimerOutput()
    sol = timed_groundstate(to, model, discretization, alg; maxiter = 100, name = "Sodium")
    return sol, to, discretization.Nₕ
end

open(joinpath(@__DIR__, "timeroutputs_profile_output.txt"), "w") do io
    for Nmesh in (20, 60)
        sol, to, Nh = profile_case(Nmesh)

        println(io, "="^90)
        println(io, "CASE: Nmesh = $Nmesh  (Nₕ = $Nh, niter = $(sol.niter), lh = 2, nspin = 1)")
        println(io, "System: Sodium, Slater exchange (spin-unpolarized)")
        println(io, "="^90)
        println(io)
        println(io, "--- Hierarchical breakdown (declaration order = structure of one SCF step) ---")
        show(io, to; compact = false, allocations = true)
        println(io)
        println(io)
        println(io, "--- Same data, sorted by total time within each level (bottleneck finder) ---")
        show(io, to; compact = false, allocations = true, sortby = :time)
        println(io)
        println(io)
    end
end

println("Wrote benchmark/timeroutputs_profile_output.txt")
