include("../../benchmarktools/atoms/setup.jl")
using DoubleFloats


using AtomicKohnSham

# LOG CONFIG
logconfig = LogConfig(orbitals_energy = true, occupation_number = true, energy = true)


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

# RESOLUTION
sol = groundstate(problem)

plot_stopping_criteria([sol]) 