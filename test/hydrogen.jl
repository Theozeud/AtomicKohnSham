@testset "Hydrogen RHF Spin Nopoloarized" begin

    logconfig = LogConfig(;orbitals_energy = false,
                           orbitals = false,
                           density = false,
                           occupation_number = false)
    problem = AtomProblem(;
                           T = Float64,
                           z = 1,
                           N = 1,
                           hartree = 1,
                           ex = NoFunctional(1),
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

    @test abs(sol.energies[:Ekin] - 0.243964867158035) < 1e-10
    @test abs(sol.energies[:Ecou] + 0.6856721558856405) < 1e-10
    @test abs(sol.energies[:Ehar] - 0.19774242154218047) < 1e-10
    @test abs(sol.energies[:Etot] + 0.24396486718542498) < 1e-10

    occ = sol.datas.occupation_number[1]
    @test occ[1] == "1s"
    @test abs(occ[2] + 0.046222445699968535) < 1e-10
    @test occ[3] == 1.0
end

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
    @test occ[3] == 0.5

    occ = sol.datas.occupation_number[2]
    @test occ[1] == "1s↓"
    @test abs(occ[2] + 0.046222445699968535) < 1e-10
    @test occ[3] == 0.5
end

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

@testset "Hydrogen Slater Spin Polarized" begin

    logconfig = LogConfig(;orbitals_energy = false,
                           orbitals = false,
                           density = false,
                           occupation_number = false)
    problem = AtomProblem(;
                           T = Float64,
                           z = 1,
                           N = 1,
                           hartree = 1,
                           ex = Functional(:lda_x, n_spin = 2),
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

    @test abs(sol.energies[:Ekin] - 0.4065340797823483) < 1e-9
    @test abs(sol.energies[:Ecou] + 0.9000752049269123) < 1e-9
    @test abs(sol.energies[:Ehar] - 0.2749225031615249) < 1e-9
    @test abs(sol.energies[:Eexc] + 0.18791545723754685) < 1e-9
    @test abs(sol.energies[:Etot] + 0.40653407922058604 ) < 1e-9

    occ = sol.datas.occupation_number[1]
    @test occ[1] == "1s↑"
    @test abs(occ[2] + 0.1942500621181194) < 1e-9
    @test occ[3] == 0.5

    occ = sol.datas.occupation_number[2]
    @test occ[1] == "1s↓"
    @test abs(occ[2] + 0.1942500621181194) < 1e-9
    @test occ[3] == 0.5
end
