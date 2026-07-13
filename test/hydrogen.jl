function run_solve(ex::TEX, ec::TEC,
                   Rmax::Real, Nmesh::Int, ordermax::Int,
                   maxiter::Int, scftol::Real) where {TEX, TEC}
    # Creation of the model
    model = KSEModel(; Z = 1, N = 1, ex = ex, ec = ec)

    # Creation of the mesh
    mesh = expmesh(0, Rmax, Nmesh; s = 1.2)

    # Creation of the Basis
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)

    # Creation of the discretization
    dis = KSEDiscretization(basis, model; lh = 0, nh=1)

    # Creation of the algorithm
    aufbau = OptimizedAufbau(max_degen = 1)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = scftol)

    # Resolution of the problem
    sol = groundstate(model, dis, alg; maxiter = maxiter)
    return sol
end

# ===================================================================
#                             SANITY CHECKS
# ===================================================================
"""
    check_sanity(sol; tol=1e-6, virial_tol=nothing)

Assert the Tier-1 physical consistency checks in `sol.sanity` (energy
balance, electron count, density normalization, orbital orthonormality) --
identities/conservation laws that must hold regardless of the system or
functional, so a violation means a real computational problem even though the
SCF loop reported `SUCCESS`.

If `virial_tol` is given, also assert the quantum virial theorem
`Ekin ≈ -Etot` (Tier 2b). This only holds when there is no exchange–
correlation functional (Hartree and nuclear-Coulomb are both homogeneous of
degree -1 in `r`, XC generally is not), and only up to the domain-truncation
defect that vanishes as `Rmax → ∞` -- so it must never be asserted for a
functional-bearing model, and needs a looser tolerance than the Tier-1 checks.
"""
function check_sanity(sol; tol::Real = 1e-6, virial_tol::Union{Real, Nothing} = nothing)
    sc = sol.sanity
    @test abs(sc.energy_balance) < tol
    @test abs(sc.electron_count_error) < tol
    @test abs(sc.density_norm_error) < tol
    @test abs(sc.orthonormality_error) < tol
    if !isnothing(virial_tol)
        @test abs(sc.virial_ratio - 1) < virial_tol
    end
end

@testset "Hydrogen RHF Spin Nopoloarized" begin

    # MODEL
    ex = NoFunctional(n_spin=1)
    ec = NoFunctional(n_spin=1)
    # DISCRETIZATION : P1IntLegendreBasis + expmesh
    Rmax        = 1000
    Nmesh       = 30
    ordermax    = 10
    # ALGORITHM : ODA + Optimized Aufbau
    maxiter     = 100
    scftol      = 1e-10

    @test_nowarn run_solve(ex, ec, Rmax, Nmesh, ordermax, maxiter, scftol)
    sol = run_solve(ex, ec, Rmax, Nmesh, ordermax, maxiter, scftol)

    @test abs(sol.energies.Ekin - 0.2439648717067077) < 1e-7
    @test abs(sol.energies.Ecou + 0.6856721626901766) < 1e-7
    @test abs(sol.energies.Ehar - 0.19774242379804155) < 1e-7
    @test abs(sol.energies.Etot + 0.24396486718542498) < 1e-7

    occ = sol.occupied[1]
    @test occ[1] == "1s"
    @test abs(occ[2] + 0.046222) < 1e-5
    @test occ[3] == 1.0

    # No exchange-correlation here (only nuclear Coulomb + Hartree, both ~1/r):
    # the virial theorem Ekin = -Etot applies, up to the Rmax-truncation defect.
    check_sanity(sol; virial_tol = 1e-5)
end


@testset "Hydrogen Slater Spin Nopolarized" begin

    # MODEL
    ex = Functional(:lda_x, n_spin=1)
    ec = NoFunctional(n_spin=1)
    # DISCRETIZATION : P1IntLegendreBasis + expmesh
    Rmax        = 1000
    Nmesh       = 30
    ordermax    = 10
    # ALGORITHM : ODA + Optimized Aufbau
    maxiter     = 100
    scftol      = 1e-10

    @test_nowarn run_solve(ex, ec, Rmax, Nmesh, ordermax, maxiter, scftol)
    sol = run_solve(ex, ec, Rmax, Nmesh, ordermax, maxiter, scftol)

    @test abs(sol.energies.Ekin - 0.4065340797823483) < 1e-7
    @test abs(sol.energies.Ecou + 0.9000752049269123) < 1e-7
    @test abs(sol.energies.Ehar - 0.2749225031615249) < 1e-7
    @test abs(sol.energies.Eexc + 0.18791545723754685) < 1e-7
    @test abs(sol.energies.Etot + 0.40653407922058604 ) < 1e-7

    occ = sol.occupied[1]
    @test occ[1] == "1s"
    @test abs(occ[2] + 0.194250) < 1e-6
    @test occ[3] == 1.0

    # Exchange functional is present here, so no virial identity to check.
    check_sanity(sol)
end


@testset "Hydrogen Perdew Spin Polarized" begin

    # MODEL
    ex = Functional(:lda_x, n_spin=2)
    ec = Functional(:lda_c_pw, n_spin=2)
    # DISCRETIZATION : P1IntLegendreBasis + expmesh
    Rmax        = 1000
    Nmesh       = 30
    ordermax    = 10
    # ALGORITHM : ODA + Optimized Aufbau
    maxiter     = 100
    scftol      = 1e-10

    @test_nowarn run_solve(ex, ec, Rmax, Nmesh, ordermax, maxiter, scftol)
    sol = run_solve(ex, ec, Rmax, Nmesh, ordermax, maxiter, scftol)

    check_sanity(sol)
end

@testset "Hydrogen Bare Coulomb (exact)" begin
    # No Hartree (hartree=0) and no exchange-correlation: the bare
    # one-electron Coulomb problem, which is exactly solvable in closed form.
    # Unlike every other testset here, this is not a regression check against
    # a value this code happened to produce in the past -- it's checked
    # against E_nl = -Z^2/(2n^2), a fact about hydrogen that holds
    # independently of this codebase.
    ex = NoFunctional(n_spin = 1)
    ec = NoFunctional(n_spin = 1)
    model = KSEModel(; Z = 1, N = 1, hartree = 0, ex = ex, ec = ec)
    mesh = expmesh(0, 1000, 30; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = 10)
    dis = KSEDiscretization(basis, model; lh = 1, nh = 2)
    aufbau = OptimizedAufbau(max_degen = 1)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-12)

    @test_nowarn groundstate(model, dis, alg; maxiter = 200)
    sol = groundstate(model, dis, alg; maxiter = 200)
    @test sol.success == "SUCCESS"

    # E_1s = -Z^2/2 exactly.
    @test isapprox(sol.energies.Etot, -0.5; atol = 1e-8)

    # sol.ϵ[l+1, k] = E_{n=k+l, l}; 2s and 2p must match the accidental
    # l-degeneracy of the pure Coulomb problem exactly.
    @test isapprox(sol.ϵ[1, 1], -1/2;  atol = 1e-8)   # 1s (n=1)
    @test isapprox(sol.ϵ[1, 2], -1/8;  atol = 1e-8)   # 2s (n=2)
    @test isapprox(sol.ϵ[2, 1], -1/8;  atol = 1e-8)   # 2p (n=2)
    @test isapprox(sol.ϵ[2, 2], -1/18; atol = 1e-8)   # 3p (n=3)

    # Pure Coulomb: the virial theorem Ekin = -Etot must hold almost exactly
    # (Rmax=1000 makes the domain-truncation defect negligible here).
    check_sanity(sol; tol = 1e-8, virial_tol = 1e-8)
end
