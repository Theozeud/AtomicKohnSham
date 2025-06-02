include("../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig(energy = false)

# CONVERGENCE WITH RESPECT TO NMESH

problem = AtomProblem(;
                T               = Float64, 
                lh              = 3, 
                method          = CDA(0.2), 
                model           = ReducedHartreeFock(54, 54), 
                Rmax            = 150.0, 
                Nmesh           = 100,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "test", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)
                        
error = convergenceNmesh(2 .^(2:7), [problem])

convergencePlotNmesh(error)



