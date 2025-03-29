include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0,
                nh              = 1, 
                method          = ODA(0.4), 
                model           = SlaterXα(1, 1), 
                Rmax            = 100, 
                Nmesh           = 10,
                typemesh        = expmesh, 
                optsmesh        = (s = 2.0,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 10,), 
                typediscre      = LDADiscretization,
                name            = "test", 
                scftol          = 1e-10,
                maxiter         = 50,
                hartree         = true,
                degen_tol       = 1e-3,
                logconfig       = logconfig,
                verbose         = 3)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])