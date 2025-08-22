include("../../benchmarktools/atoms/setup.jl")

using DoubleFloats
using AtomicKohnSham
using Libxc

# LOG CONFIG
logconfig = LogConfig(; orbitals_energy = true,
                        orbitals = false,
                        density = false,
                        occupation_number = true
                        )


problem = AtomProblem(;
                        T =Double64,
                        z = 26,
                        N = 26,
                        hartree = 1,
                        ex = NoFunctional(1),#Functional(:lda_x, n_spin = 2),
                        ec = NoFunctional(1),#Functional(:lda_c_pw, n_spin = 2),
                        lh = 2,
                        nh = 5,
                        Nmesh =  20,
                        Rmax = 10000,
                        typemesh = explinmesh,
                        optsmesh = (s = 2.0,rswitch=800, nlin=10),
                        typebasis = P1IntLegendreBasis,
                        optsbasis = (ordermax = 14,),
                        integration_method = GaussLegendre,
                        optsintegration = (npoints = 2000,),
                        alg = ODA(0.4),
                        scftol = 1e-15,
                        maxiter = 200,
                        degen_tol = 1e-2,
                        logconfig = logconfig,
                        verbose = 1)

# RESOLUTION
@time sol = groundstate(problem);
print(sol)

# PLOT STOPPING CRITERIA
p = plot_stopping_criteria([sol])

# PLOT DENSITY
X = AtomicKohnSham.exprange(0.001,9999.5, 10000; s = 1.2)
ρ = plot_density([sol], X)

# PLOT EIGENVECTOR
E = plot_eigenvector(sol,X)
