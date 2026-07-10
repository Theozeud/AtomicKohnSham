# Compares the eigensolver currently used in `find_orbital!`
# (`AtomicKohnSham/src/algorithms/oda/loop.jl`) — a full dense
# `eigen(Symmetric(A))` that discards all but the lowest `nₕ` eigenpairs —
# against a range-restricted `eigen(Symmetric(A), 1:nₕ)`, which asks LAPACK
# directly for only the eigenpairs actually used downstream.
#
# The comparison runs on the *real* Kohn-Sham Hamiltonian matrices produced
# by converged SCF solves (not random matrices: a range-restricted
# eigensolver's performance depends on the eigenvalue distribution, and a
# physical Hamiltonian's shell-clustered spectrum behaves very differently
# from a random Wigner matrix).
#
# Run with: julia benchmark/eigensolver_comparison.jl

cd(@__DIR__)
using Pkg
Pkg.activate(".")

using AtomicKohnSham
using Libxc
using LinearAlgebra
using Statistics
using BenchmarkTools
using Printf

"""
    converged_discretization(Nmesh; Z, N, Rmax, ordermax, lh) -> KSEDiscretization

Solve a spin-unpolarized Slater-exchange ground state and return its
discretization, whose `ksham.H` blocks hold the final (converged) Kohn-Sham
Hamiltonian — the same matrices `find_orbital!` diagonalizes every SCF
iteration.
"""
function converged_discretization(Nmesh::Int; Z::Real = 11, N::Real = 11,
                                 Rmax::Real = 500, ordermax::Int = 10, lh::Int = 2)
    model = KSEModel(; Z = Z, N = N, ex = Functional(:lda_x; n_spin = 1), ec = NoFunctional(1))
    mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
    dis = KSEDiscretization(basis, model; lh = lh, nh = 10,
        fem_integration_method = GaussLegendre(basis, 2000))
    aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-1)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-9)
    groundstate(model, dis, alg; maxiter = 100)
    return dis
end

"Reconstruct the exact matrix `find_orbital!` diagonalizes for channel `(l, σ)`."
function hamiltonian_matrix(dis::KSEDiscretization, l::Int; σ::Int = 1)
    S = dis.femops.S
    @views vH = dis.ksham.H[:, :, l + 1, σ]
    return Symmetric(S * vH * S)
end

"""
    compare_channel(A, nh) -> (t_full_ms, t_part_ms, speedup, max_abs_error)

Benchmark full vs. range-restricted `eigen` on Hamiltonian block `A`, and
check that the lowest `nh` eigenvalues agree between the two methods.
"""
function compare_channel(A::Symmetric, nh::Int)
    b_full = @benchmark eigen($A) samples=30 evals=1
    b_part = @benchmark eigen($A, 1:$nh) samples=30 evals=1

    λ_full = eigen(A).values[1:nh]
    λ_part = eigen(A, 1:nh).values
    err = maximum(abs.(λ_full .- λ_part))

    t_full = median(b_full).time / 1e6   # ns -> ms
    t_part = median(b_part).time / 1e6
    return t_full, t_part, t_full / t_part, err
end

function run_comparison(; Nmeshes = (20, 40, 60), lh::Int = 2)
    println("Comparing eigen(Symmetric(A)) vs eigen(Symmetric(A), 1:nh)")
    println("on real converged Kohn-Sham Hamiltonians (Sodium, Slater exchange)\n")

    for Nmesh in Nmeshes
        dis = converged_discretization(Nmesh; lh = lh)
        Nh, nh = dis.Nₕ, dis.nₕ

        println("="^78)
        @printf("Nmesh = %d  ->  Nh = %d, nh = %d\n", Nmesh, Nh, nh)
        println("="^78)
        @printf("%-6s %12s %12s %10s %12s\n", "l", "full (ms)", "partial (ms)", "speedup", "max |Δε|")
        for l in 0:lh
            A = hamiltonian_matrix(dis, l)
            t_full, t_part, speedup, err = compare_channel(A, nh)
            @printf("%-6d %12.4f %12.4f %9.2fx %12.2e\n", l, t_full, t_part, speedup, err)
        end
        println()
    end
end

run_comparison()
