include("../../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# MODEL
z = 8
N = 8

# LOG CONFIG
logconfig = LogConfig()

problem = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(z, N), 
                Rmax            = 80.0, 
                Nmesh           = 100,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = ShortP1IntLegendreBasis, 
                optsbasis       = (ordermax = 5,), 
                name            = "test", 
                scftol          = 1e-8,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)


# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol])