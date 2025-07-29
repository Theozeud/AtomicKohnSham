include("../../benchmarktools/atoms/setup.jl")
using DoubleFloats


using AtomicKohnSham
using Libxc

# LOG CONFIG
logconfig = LogConfig(; orbitals_energy = false, 
                        orbitals = false, 
                        density = false, 
                        occupation_number = false)

model = KSEModel(;z=15,N=15,ex=Functional(:lda_x, n_spin=2), ec=Functional(:lda_c_pw, n_spin=2))

#model = KSEModel(;z=26,N=26,ex=Functional(:lda_x, n_spin=1), ec=NoFunctional(1))


problem = AtomProblem(;
                T               = Float64,
                lh              = 1,
                nh              = 3,
                alg             = ODA(0.4),
                model           = model,
                Rmax            = 500,
                Nmesh           = 50,
                typemesh        = expmesh,
                optsmesh        = (s = 1.5,),
                typebasis       = P1IntLegendreBasis,
                optsbasis       = (ordermax = 20,),
                name            = "Hydrogen",
                scftol          = 1e-11,
                maxiter         = 50,
                degen_tol       = 1e-2,
                logconfig       = logconfig,
                verbose         = 3)

# RESOLUTION
@time sol = groundstate(problem);

plot_stopping_criteria([sol])


# TESTS

using AtomicKohnSham: init_cache!
mesh = expmesh(0, 200.0, 40; s=1.5)
basis = P1IntLegendreBasis(mesh; ordermax=10)

d = KSEDiscretization(1, basis, mesh, 1, 2, GaussLegendre(basis,100))
init_cache!(d,model)

plot_density([sol], exprange(0,100,1000;s=1.5))