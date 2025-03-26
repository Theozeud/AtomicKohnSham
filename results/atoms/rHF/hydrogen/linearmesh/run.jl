include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

problem = AtomProblem(;
                T               = Float64, 
                lh              = 1,
                nh              = 100, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(8, 8), 
                Rmax            = 80, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s = 0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                typediscre      = LDADiscretization,
                name            = "test", 
                scftol          = 1e-10,
                maxiter         = 20,
                hartree         = true,
                degen_tol       = 1e-10,
                logconfig       = logconfig,
                verbose         = 1)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])