module AtomicKohnSham

    using LinearAlgebra
    using SparseArrays
    using TensorOperations
    using FastGaussQuadrature
    using QuadGK: gauss as quadgk_gauss
    using KrylovKit
    using Arpack
    using Optim
    using Integrals
    using HypergeometricFunctions
    using UnPack
    using Base.Threads

    using Libxc: Functional, OptArray
    import Libxc: is_lda
    import Libxc: evaluate as evaluate_functional
    import Libxc: evaluate! as evaluate_functional!


    # =====================================
    #               DATA
    # =====================================
    export ATOMIC_NUMBER_TO_NAME, ATOMIC_NUMBER_TO_SYMBOL
    include("data/periodic_tables.jl")

    # =====================================
    #               UTILS
    # =====================================
    include("utils/free allocation.jl")
    include("utils/sparse tools.jl")
    include("utils/orbital_labels.jl")

    # =====================================
    #               POLYSET
    # =====================================
    # ALL THESE TOOLS SHOULD MOVE TO A NEW PACKAGE POLYSET
    export PolySet, Legendre, IntLegendre
    export evaluate, evaluate!, integrate!, allocate_polyset, mul, npolys, maxdeg

    include("fem/polynomial.jl")
    include("fem/legendre polynomial.jl")

    # =====================================
    #               FEM TOOLS
    # =====================================
    export FEMBasis

    export P1IntLegendreGenerator, P1IntLegendreBasis

    export mass_matrix, sparse_mass_matrix,
        stiffness_matrix, sparse_stiffness_matrix,
        mass_tensor, mass_vector

    export FEMIntegrationMethod, ExactIntegration, QuadratureIntegration, GaussLegendre

    abstract type AbstractGenerator{T} end
    @inline Base.eltype(::AbstractGenerator{T}) where {T} = T
    @inline Base.length(gen::AbstractGenerator) = gen.size
    @inline getpolynomial(gen::AbstractGenerator, n::Int) = gen[n]
    @inline getderivpolynomial(gen::AbstractGenerator, n::Int) = getderivpolynomial(gen)[n]

    export Mesh, linmesh, geometricmesh, polynomialmesh, expmesh, explinmesh
    include("fem/mesh.jl")

    include("fem/generators.jl")
    include("fem/basis.jl")
    include("fem/integration methods.jl")
    include("fem/weights.jl")
    include("fem/matrices.jl")
    include("fem/local matrix.jl")
    include("fem/integration_formula.jl")

    # =====================================
    #              PHYSICS
    # =====================================
    export NoFunctional, BuiltinFunctional
    export evaluate_functional, evaluate_functional!
    export KSEModel, RHF, Slater

    include("physics/exchange correlation/BuiltinFunctional.jl")
    include("physics/exchange correlation/SlaterXa.jl")
    include("physics/models.jl")

    # =====================================
    #               DISCRETIZATION
    # =====================================
    export KSEDiscretization
    include("discretization/discretization.jl")
    include("discretization/assemble.jl")
    include("discretization/energies.jl")
    include("discretization/density.jl")

    abstract type SCFAlgorithm end
    abstract type SCFCache end

    # =====================================
    #               LogBook
    # =====================================
    export LogBook
    include("solver/logbook.jl")


    export Energies
    include("solver/energies.jl")

    # =====================================
    #               ALGORITHMS
    # =====================================
    export line_search_energy
    include("algorithms/optimization/line_search.jl")

    abstract type Aufbau end
    abstract type AufbauCache end
    export OptimizedAufbau, SmearedAufbau, FrozenAufbau
    const AUFBAU_METHOD = [:optimized, :smeared, :frozen]
    export aufbau!
    include("algorithms/aufbau/optimized.jl")
    include("algorithms/aufbau/smeared.jl")
    include("algorithms/aufbau/frozen.jl")


    export KSESolver
    include("solver/solver.jl")

    export CallbackSet
    include("solver/callback.jl")


    export ODA
    include("algorithms/oda/types.jl")
    include("algorithms/oda/loop.jl")


    # =====================================
    #               SOLVER
    # =====================================
    export KSEContext
    include("solution/context.jl")

    export KSESolution
    export eval_orbital, eval_density, eval_hartree, eval_nuclear,
        eval_kinetic_potential, eval_vxc, eval_effective_potential
    export SanityChecks
    include("solution/sanity.jl")
    include("solution/types.jl")
    include("solution/show.jl")
    include("solution/postprocess.jl")

    export write_report
    include("solution/report.jl")

    export groundstate
    include("solver/groundstate.jl")

    export LogFileCallback, write_log_header
    include("solver/logging.jl")

    # =====================================
    #             PLOTTING (extension)
    # =====================================
    export plot_density, plot_orbitals, plot_potentials, plot_convergence,
        plot_energy_breakdown
    include("solution/plotting_stubs.jl")
end
