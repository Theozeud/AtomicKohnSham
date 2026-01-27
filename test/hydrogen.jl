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

    @test abs(sol.energies.Ekin - 0.2439648717067077) < 1e-10
    @test abs(sol.energies.Ecou + 0.6856721626901766) < 1e-10
    @test abs(sol.energies.Ehar - 0.19774242379804155) < 1e-10
    @test abs(sol.energies.Etot + 0.24396486718542498) < 1e-10

    occ = sol.occupied[1]
    @test occ[1] == "1s"
    @test abs(occ[2] + 0.046222443387387145) < 1e-10
    @test occ[3] == 1.0
end

#=
@testset "Hydrogen RHF Spin Polarized" begin

    logconfig = LogConfig(;orbitals_energy = false,
                           orbitals = false,
                           density = false,
                           occupation_number = false)
    problem = AtomProblem(;
                           T = Float64,
                           z = 1,
                           N = 1,
                           hartree = 1,
                           ex = NoFunctional(2),
                           ec = NoFunctional(2),
                           lh = 0,
                           nh = 1,
                           Nmesh = 10,
                           Rmax = 300,
                           typemesh = expmesh,
                           optsmesh = (s = 1.5,),
                           typebasis = P1IntLegendreBasis,
                           optsbasis = (ordermax = 10,),
                           integration_method = GaussLegendre,
                           optsintegration = (npoints = 1000,),
                           alg = ODA(0.4),
                           scftol = 1e-11,
                           maxiter = 100,
                           degen_tol = 1e-2,
                           logconfig = logconfig,
                           verbose = 0)



    @test_nowarn sol = groundstate(problem)
    sol = groundstate(problem)

    @test abs(sol.energies[:Ekin] - 0.243964867158035) < 1e-10
    @test abs(sol.energies[:Ecou] + 0.6856721558856405) < 1e-10
    @test abs(sol.energies[:Ehar] - 0.19774242154218047) < 1e-10
    @test abs(sol.energies[:Etot] + 0.24396486718542498) < 1e-10

    occ = sol.datas.occupation_number[1]
    @test occ[1] == "1s↑"
    @test abs(occ[2] + 0.046222445699968535) < 1e-10
    #@test occ[3] == 0.5

    occ = sol.datas.occupation_number[2]
    @test occ[1] == "1s↓"
    @test abs(occ[2] + 0.046222445699968535) < 1e-10
    #@test occ[3] == 0.5
end
=#

#=
@testset "Hydrogen Slater Spin Nopolarized" begin

    logconfig = LogConfig(;orbitals_energy = false,
                           orbitals = false,
                           density = false,
                           occupation_number = false)
    problem = AtomProblem(;
                           T = Float64,
                           z = 1,
                           N = 1,
                           hartree = 1,
                           ex = Functional(:lda_x, n_spin = 1),
                           ec = NoFunctional(1),
                           lh = 0,
                           nh = 1,
                           Nmesh = 10,
                           Rmax = 300,
                           typemesh = expmesh,
                           optsmesh = (s = 1.5,),
                           typebasis = P1IntLegendreBasis,
                           optsbasis = (ordermax = 10,),
                           integration_method = GaussLegendre,
                           optsintegration = (npoints = 1000,),
                           alg = ODA(0.4),
                           scftol = 1e-11,
                           maxiter = 100,
                           degen_tol = 1e-2,
                           logconfig = logconfig,
                           verbose = 0)



    @test_nowarn sol = groundstate(problem)
    sol = groundstate(problem)

    @test abs(sol.energies[:Ekin] - 0.4065340797823483) < 1e-9
    @test abs(sol.energies[:Ecou] + 0.9000752049269123) < 1e-9
    @test abs(sol.energies[:Ehar] - 0.2749225031615249) < 1e-9
    @test abs(sol.energies[:Eexc] + 0.18791545723754685) < 1e-9
    @test abs(sol.energies[:Etot] + 0.40653407922058604 ) < 1e-9

    occ = sol.datas.occupation_number[1]
    @test occ[1] == "1s"
    @test abs(occ[2] + 0.1942500621181194) < 1e-9
    @test occ[3] == 1.0
end
=#

#=
@testset "Hydrogen Slater Spin Polarized" begin


    Model = ModelParams()

    problem = AtomProblem(;
        # PRECISION
        T = Float64,
        # MODEL
        z = 1,
        N = 1,
        hartree = 1,
        ex = Functional(:lda_x, n_spin = 2),
        ec = NoFunctional(2),

        # DISCRETIZATION
        lh = 0,
        nh = 1,
        Nmesh = 10,
        Rmax = 300,
        mesh = :exp,
        optsmesh = (s = 1.5,),
        basis = :p1intlegendre,
        optsbasis = (ordermax = 10,),
        integration_method = GaussLegendre,
        optsintegration = (npoints = 1000,),

        # ALGORITHME
        alg = ODA(0.4),
        scftol = 1e-11,
        maxiter = 100,
        degen_tol = 1e-2,
        verbose = 0)



    @test_nowarn sol = groundstate(problem)
    sol = groundstate(problem)

    @test abs(sol.energies[:Ekin] - 0.4065340797823483) < 1e-9
    @test abs(sol.energies[:Ecou] + 0.9000752049269123) < 1e-9
    @test abs(sol.energies[:Ehar] - 0.2749225031615249) < 1e-9
    @test abs(sol.energies[:Eexc] + 0.18791545723754685) < 1e-9
    @test abs(sol.energies[:Etot] + 0.40653407922058604 ) < 1e-9

    occ = sol.datas.occupation_number[1]
    @test occ[1] == "1s↑"
    @test abs(occ[2] + 0.1942500621181194) < 1e-9
    #@test occ[3] == 0.5

    occ = sol.datas.occupation_number[2]
    @test occ[1] == "1s↓"
    @test abs(occ[2] + 0.1942500621181194) < 1e-9
    #@test occ[3] == 0.5
end
=#
