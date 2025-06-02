include("../../../../../benchmarktools/atoms/setup.jl")
using AtomicKohnSham

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)

using DoubleFloats

problem = AtomProblem(;
                T               = Float64, 
                lh              = 3,
                nh              = 5, 
                alg             = ODA(0.4), 
                model           = RHF(;z=30, N=30),  
                Rmax            = 1000, 
                Nmesh           = 50,
                typemesh        = expmesh, 
                optsmesh        = (s = 1.5,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 20,), 
                name            = "test", 
                scftol          = 1e-13,
                maxiter         = 60,
                degen_tol       = 1e-2,
                logconfig       = logconfig,
                verbose         = 1)

# RESOLUTION
@time sol = groundstate(problem)

plot_stopping_criteria([sol]) 