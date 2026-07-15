# Stage 4 (see PLAN.md): end-to-end SCF smoke test for GGA (PBE) on Hydrogen
# and Helium, mirroring test/hydrogen.jl's run_solve/check_sanity pattern.
# Only Tier-1 sanity checks are asserted (energy balance, electron count,
# density norm, orthonormality) -- not the virial ratio, since GGA XC (like
# LDA XC) breaks the pure-Coulomb virial identity (see check_sanity's
# docstring in test/hydrogen.jl).
#
# Run: julia --project=. dev/gga/stage4_endtoend.jl

using AtomicKohnSham
using Test

function run_solve_gga(Z::Real, N::Real, nspin::Int,
                       Rmax::Real, Nmesh::Int, ordermax::Int,
                       maxiter::Int, scftol::Real)
    model = PBE(Z = Z, N = N, nspin = nspin)
    mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
    dis = KSEDiscretization(basis, model; lh = 0, nh = 1)
    aufbau = OptimizedAufbau(max_degen = 1)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = scftol)
    return groundstate(model, dis, alg; maxiter = maxiter)
end

function check_sanity(sol; tol::Real = 1e-6)
    sc = sol.sanity
    @test abs(sc.energy_balance) < tol
    @test abs(sc.electron_count_error) < tol
    @test abs(sc.density_norm_error) < tol
    @test abs(sc.orthonormality_error) < tol
end

@testset "Hydrogen PBE (nspin=1)" begin
    sol = run_solve_gga(1.0, 1.0, 1, 1000, 30, 10, 100, 1e-10)
    println("Hydrogen PBE: success=$(sol.success)  Etot=$(sol.energies.Etot)")
    @test sol.success == "SUCCESS"
    check_sanity(sol)
end

@testset "Helium PBE (nspin=1)" begin
    sol = run_solve_gga(2.0, 2.0, 1, 1000, 30, 10, 100, 1e-10)
    println("Helium PBE: success=$(sol.success)  Etot=$(sol.energies.Etot)")
    @test sol.success == "SUCCESS"
    check_sanity(sol)
end
