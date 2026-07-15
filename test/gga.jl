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

# check_sanity is defined in test/hydrogen.jl (included before this file in
# runtests.jl); GGA XC breaks the pure-Coulomb virial identity just like LDA
# does, so virial_tol is never passed here.

@testset "Hydrogen PBE Spin Nopolarized" begin
    Rmax        = 1000
    Nmesh       = 30
    ordermax    = 10
    maxiter     = 100
    scftol      = 1e-10

    @test_nowarn run_solve_gga(1.0, 1.0, 1, Rmax, Nmesh, ordermax, maxiter, scftol)
    sol = run_solve_gga(1.0, 1.0, 1, Rmax, Nmesh, ordermax, maxiter, scftol)

    @test abs(sol.energies.Ekin - 0.44724800200785514) < 1e-6
    @test abs(sol.energies.Ecou + 0.9446278659341152) < 1e-6
    @test abs(sol.energies.Ehar - 0.28981159147814184) < 1e-6
    @test abs(sol.energies.Eexc + 0.2513603248087004) < 1e-6
    @test abs(sol.energies.Etot + 0.45892859725681856) < 1e-6

    occ = sol.occupied[1]
    @test occ[1] == "1s"
    @test occ[3] == 1.0

    check_sanity(sol)
end

@testset "Helium PBE Spin Nopolarized" begin
    Rmax        = 1000
    Nmesh       = 30
    ordermax    = 10
    maxiter     = 100
    scftol      = 1e-10

    @test_nowarn run_solve_gga(2.0, 2.0, 1, Rmax, Nmesh, ordermax, maxiter, scftol)
    sol = run_solve_gga(2.0, 2.0, 1, Rmax, Nmesh, ordermax, maxiter, scftol)

    @test abs(sol.energies.Ekin - 2.8559475957199343) < 1e-6
    @test abs(sol.energies.Ecou + 6.729457535623727) < 1e-6
    @test abs(sol.energies.Ehar - 2.0267373340851425) < 1e-6
    @test abs(sol.energies.Eexc + 1.046162261029548) < 1e-6
    @test abs(sol.energies.Etot + 2.892934866848198) < 1e-6

    occ = sol.occupied[1]
    @test occ[1] == "1s"
    @test occ[3] == 2.0

    check_sanity(sol)
end

@testset "Mixed-rung exchange-correlation is rejected" begin
    ex = Functional(:gga_x_pbe, n_spin = 1)
    ec = BuiltinFunctional(:lda_c_pw; nspin = 1)
    @test_throws ErrorException KSEModel(; Z = 1.0, N = 1.0, ex = ex, ec = ec)
end
