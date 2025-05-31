include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

using DoubleFloats

problem = AtomProblem(;
                T               = Float64, 
                lh              = 1,
                nh              = 2, 
                method          = ODA(0.4), 
                model           = SlaterXα(10, 10),  
                Rmax            = 100, 
                Nmesh           = 10,
                typemesh        = expmesh, 
                optsmesh        = (s = 1.5,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 20,), 
                typediscre      = LDADiscretization,
                name            = "test", 
                scftol          = 1e-13,
                maxiter         = 60,
                hartree         = true,
                degen_tol       = 1e-2,
                logconfig       = logconfig,
                verbose         = 1)

# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol]) 