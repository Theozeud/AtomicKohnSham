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
    sol = groundstate(model, dis, alg; maxiter = maxiter, verbose = 0)
    return sol
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
end
