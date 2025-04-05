include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0,
                nh              = 1, 
                method          = ODA(0.4), 
                model           = LSDA(1, 1),  
                Rmax            = 300, 
                Nmesh           = 10,
                typemesh        = expmesh, 
                optsmesh        = (s = 1.5,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 15,), 
                typediscre      = LSDADiscretization,
                name            = "test", 
                scftol          = 1e-8,
                maxiter         = 10,
                hartree         = true,
                degen_tol       = 1e-2,
                logconfig       = logconfig,
                verbose         = 3)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])