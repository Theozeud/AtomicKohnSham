
@testset "Oxygen RHF" begin

    @test_nowarn logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)
    logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

    @test_nowarn problem = AtomProblem(;
                    T               = Float64, 
                    lh              = 1,
                    nh              = 2, 
                    alg             = ODA(0.4), 
                    model           = RHF(;z=8, N=8),  
                    Rmax            = 100, 
                    Nmesh           = 10,
                    typemesh        = expmesh, 
                    optsmesh        = (s = 1.5,), 
                    typebasis       = P1IntLegendreGenerator, 
                    optsbasis       = (ordermax = 20,), 
                    name            = "Oxygen", 
                    scftol          = 1e-11,
                    maxiter         = 60,
                    degen_tol       = 1e-2,
                    logconfig       = logconfig,
                    verbose         = 0)

    problem = AtomProblem(;
                    T               = Float64, 
                    lh              = 1,
                    nh              = 2, 
                    alg             = ODA(0.4), 
                    model           = RHF(;z=8, N=8),  
                    Rmax            = 100, 
                    Nmesh           = 10,
                    typemesh        = expmesh, 
                    optsmesh        = (s = 1.5,), 
                    typebasis       = P1IntLegendreGenerator, 
                    optsbasis       = (ordermax = 20,), 
                    name            = "Oxygen", 
                    scftol          = 1e-11,
                    maxiter         = 60,
                    degen_tol       = 1e-2,
                    logconfig       = logconfig,
                    verbose         = 0)

    @test_nowarn sol = groundstate(problem)
    sol = groundstate(problem)

    @test abs(sol.energies[:Ekin] - 67.03856447670336) < 1e-9
    @test abs(sol.energies[:Ecou] + 166.0557278383923) < 1e-9
    @test abs(sol.energies[:Ehar] - 31.978598417860578) < 1e-9 
    @test abs(sol.energies[:Etot] + 67.03856494382838) < 1e-9 

    occ1 = sol.datas.occupation_number[1]
    @test occ1[1] == "1s"
    @test abs(occ1[2] + 16.91253843847191) < 1e-9
    @test occ1[3] == 2.0

    occ2 = sol.datas.occupation_number[2]
    @test occ2[1] == "2s"
    @test abs(occ2[2] + 0.5228830530863569) < 1e-9
    @test occ2[3] == 2.0

    occ3 = sol.datas.occupation_number[3]
    @test occ3[1] == "2p"
    @test abs(occ3[2] + 0.04728084046941642) < 1e-9
    @test occ3[3] == 4.0
end