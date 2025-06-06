module AtomicKohnSham

    # -------------------------------------
    #               DEPENDANCIES
    # -------------------------------------
    using LinearAlgebra
    using SparseArrays
    using FillArrays
    using BlockDiagonals

    using TensorOperations

    using FastGaussQuadrature
    using KrylovKit
    using Arpack
    using Optim
    using Integrals

    using HypergeometricFunctions

    using UnPack
    
    using Base.Threads


    
    # -------------------------------------
    #               UTILS
    # -------------------------------------
    include("utils/free allocation.jl")
    include("utils/sparse tools.jl")

    
    # -------------------------------------
    #               POLYSET
    # -------------------------------------
    # ALL THESE TOOLS SHOULD MOVE TO A NEW PACKAGE POLYSET
    export PolySet, Legendre, IntLegendre
    export evaluate!, integrate!, allocate_polyset, mul, npolys, maxdeg

    include("fem/polynomial.jl")
    include("fem/legendre polynomial.jl")


    # -------------------------------------
    #               FEM TOOLS
    # -------------------------------------
    export FEMBasis

    export P1IntLegendreGenerator

    export mass_matrix, sparse_mass_matrix, 
           stiffness_matrix, sparse_stiffness_matrix,
           mass_tensor
    
    export IntegrationMethod, ExactIntegration, QuadratureIntegration

    abstract type AbstractGenerator{T} end
    @inline Base.eltype(::AbstractGenerator{T}) where T = T
    @inline Base.length(gen::AbstractGenerator) = gen.size
    @inline getpolynomial(gen::AbstractGenerator, n::Int) = gen[n]
    @inline getderivpolynomial(gen::AbstractGenerator, n::Int) = getderivpolynomial(gen)[n]

    export Mesh, linmesh, geometricmesh, polynomialmesh, expmesh
    include("fem/mesh.jl")

    include("fem/generators.jl")
    include("fem/basis.jl")
    include("fem/computations.jl")
    include("fem/matrices.jl")
    include("fem/local matrix.jl")
    include("fem/integration_formula.jl")


    # -------------------------------------
    #               MODELISATION
    # -------------------------------------
    export ExchangeCorrelation, NoExchangeCorrelation, SlaterXα, LSDA
    export KSEModel, RHF, SlaterXα

    include("exchange correlation.jl")
    include("models.jl")
    

    # -------------------------------------
    #               DISCRETIZATION
    # -------------------------------------
    export KSEDiscretization
    include("discretization/discretization.jl")
    include("discretization/operator.jl")
    include("discretization/energies.jl")
    include("discretization/density.jl")

    
    # -------------------------------------
    #               LogBook
    # -------------------------------------
    export LogConfig, LogBook
    include("log.jl")

    # -------------------------------------
    #            SCF ALGORITHMS
    # -------------------------------------
    export KSESolver, SolverOptions
    export SCFAlgorithm, CDA, ODA, Quadratic

    include("algorithms/abstract.jl")

    include("solver.jl")

    include("algorithms/rca.jl")
    include("algorithms/oda_procedure.jl")
    #include("algorithms/quadratic.jl")
    include("algorithms/combined.jl")

    export aufbau!
    include("algorithms/aufbau.jl")

    # -------------------------------------
    #       PROBLEM-groundstate-SOLUTION
    # -------------------------------------
    export AtomProblem
    include("problem.jl")

    export KSESolution, eigenvector, density, total_charge
    include("solution.jl")

    export groundstate
    include("groundstate.jl")
end
