include("../../benchmarktools/atoms/setup.jl")

using DoubleFloats
using AtomicKohnSham
using Libxc

# LOG CONFIG
logconfig = LogConfig(; orbitals_energy = false,
                        orbitals = false,
                        density = false,
                        occupation_number = false
                        )

function density(Z::AbstractVector)

    Sols = []
    for z ∈ Z

        # CREATION OF THE PROBLEM
        lh = Int(min(ceil(z/20),3))
        nh = min(max(z,2),10)
        problem = AtomProblem(;
                                T =Float64,
                                z = z,
                                N = z,
                                hartree = 1,
                                ex = Functional(:lda_x, n_spin = 1),#NoFunctional(1), #Functional(:lda_x, n_spin = 1),
                                ec = NoFunctional(1),#Functional(:lda_c_pw, n_spin = 2),
                                lh = lh,
                                nh = nh,
                                Nmesh = 30,
                                Rmax = 1000,
                                typemesh = expmesh,
                                optsmesh = (s = 1.2,),
                                typebasis = P1IntLegendreBasis,
                                optsbasis = (ordermax = 10,),
                                integration_method = GaussLegendre,
                                optsintegration = (npoints = 1000,),
                                alg = ODA(0.4),
                                scftol = 1e-10,
                                maxiter = 100,
                                degen_tol = 1e-2,
                                logconfig = logconfig,
                                verbose = 0)

        # RESOLUTION
        @time sol = groundstate(problem);
        push!(Sols, sol)

    end

    return Sols


end

Sols = density([3, 16, 17, 18, 20, 26])

X = AtomicKohnSham.exprange(0.001,99.9, 10000; s = 1.5)
fig = plot_density(Sols, X)
save("./examples/density computations/densities.pdf",fig)
