include("../../benchmarktools/atoms/setup.jl")

using DoubleFloats
using AtomicKohnSham
using Libxc

# LOG CONFIG
logconfig = LogConfig(; orbitals_energy = false,
                        orbitals = false,
                        density = false,
                        occupation_number = false)


problem = AtomProblem(;
                        T = Float64,
                        z = 1,
                        N = 1,
                        hartree = 1,
                        ex = Functional(:lda_x, n_spin = 1),
                        ec = NoFunctional(1),#Functional(:lda_c_pw, n_spin = 2),
                        lh = 0,
                        nh = 1,
                        Nmesh = 15,
                        Rmax = 300,
                        typemesh = expmesh,
                        optsmesh = (s = 1.5,),
                        typebasis = P1IntLegendreBasis,
                        optsbasis = (ordermax = 10,),
                        integration_method = GaussLegendre,
                        optsintegration = (npoints = 1000,),
                        alg = ODA(0.4),
                        scftol = 1e-10,
                        maxiter = 50,
                        degen_tol = 1e-2,
                        logconfig = logconfig,
                        verbose = 3)

# RESOLUTION
@time sol = groundstate(problem);

plot_stopping_criteria([sol])

# TESTS
plot_density([sol], AtomicKohnSham.exprange(0, 299, 1000; s = 1.5))
