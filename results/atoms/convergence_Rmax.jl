include("../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# MODEL
z = 1
N = 1
model = ReducedHartreeFock(z, N)

# LOG CONFIG
logconfig = LogConfig(energy = false)

# CONVERGENCE WITH RESPECT TO NMESH

problem = AtomProblem(;
                T               = Float64, 
                lh              = 0, 
                method          = CDA(0.7), 
                model           = ReducedHartreeFock(z, N), 
                Rmax            = 40.0, 
                Nmesh           = 50,
                typemesh        = linmesh, 
                optsmesh        = (;), 
                typebasis       = ShortP1IntLegendreBasis, 
                optsbasis       = (ordermax = 3,), 
                name            = "test", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)
                        
error = convergenceRmax(LinRange(5,60,20), [problem])

convergencePlotRmax(error)

