
@testset "Hydrogen RHF" begin

    @test_nowarn logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)
    logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

    @test_nowarn problem = AtomProblem(;
                    T               = Float64, 
                    lh              = 0,
                    nh              = 1, 
                    alg             = ODA(0.4), 
                    model           = RHF(;z=1, N=1),  
                    Rmax            = 100, 
                    Nmesh           = 10,
                    typemesh        = expmesh, 
                    optsmesh        = (s = 1.5,), 
                    typebasis       = P1IntLegendreGenerator, 
                    optsbasis       = (ordermax = 20,), 
                    name            = "Hydrogen", 
                    scftol          = 1e-11,
                    maxiter         = 60,
                    degen_tol       = 1e-2,
                    logconfig       = logconfig,
                    verbose         = 0)

    problem = AtomProblem(;
                    T               = Float64, 
                    lh              = 0,
                    nh              = 1, 
                    alg             = ODA(0.4), 
                    model           = RHF(;z=1, N=1),  
                    Rmax            = 100, 
                    Nmesh           = 10,
                    typemesh        = expmesh, 
                    optsmesh        = (s = 1.5,), 
                    typebasis       = P1IntLegendreGenerator, 
                    optsbasis       = (ordermax = 20,), 
                    name            = "Hydrogen", 
                    scftol          = 1e-11,
                    maxiter         = 60,
                    degen_tol       = 1e-2,
                    logconfig       = logconfig,
                    verbose         = 0)

    @test_nowarn sol = groundstate(problem)
    sol = groundstate(problem)

    @test abs(sol.energies[:Ekin] - 0.24396486583278804) < 1e-10
    @test abs(sol.energies[:Ecou] + 0.685672154045772) < 1e-10
    @test abs(sol.energies[:Ehar] - 0.19774242102873807) < 1e-10 
    @test abs(sol.energies[:Etot] + 0.24396486718424593) < 1e-10 

    occ = sol.datas.occupation_number[1]
    @test occ[1] == "1s"
    @test abs(occ[2] + 0.04622244682275613) < 1e-10
    @test occ[3] == 1.0
end