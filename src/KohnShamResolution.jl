module KohnShamResolution

    # DEPENDANCIES
    using LinearAlgebra
    using SparseArrays
    using FillArrays
    using BlockDiagonals

    using TensorOperations
    using LinearOperators
    
    using FastGaussQuadrature
    using Krylov
    using Arpack
    using Optim
    using Integrals

    using HypergeometricFunctions

    using UnPack
    using Memoize
    
    using Base.Threads
    
    # ANNEXE
    include("utils.jl")
    include("maths.jl")

    # MESH
    export Mesh, linmesh, geometricmesh, polynomialmesh, expmesh
    include("mesh.jl")

    ############        LAURENT POLYNOMIAL        ############
    export AbstractPolynomial

    abstract type AbstractPolynomial{T} end
    @inline Base.eltype(::AbstractPolynomial{T}) where T = T

    export LaurentPolynomial
    export Polynomial, Monomial, RandMonomial, RandPolynomial, RootsPolynomial
    export deg, degmax, degmin, ismonomial, iszero
    export integrate, deriv, scalar_product, normL2, elag!, diveucl
    include("fem/laurentpolynomial.jl")

    export Legendre, intLegendre
    include("fem/legendre polynomial.jl")


    ############        FINITE ELEMENT METHOD        ############

    abstract type Basis end
    export PolynomialBasis

    export P1IntLegendreGenerator

    export mass_matrix, sparse_mass_matrix, 
           stiffness_matrix, sparse_stiffness_matrix,
           weight_mass_matrix, sparse_weight_mass_matrix,
           weight_mass_vector, sparse_weight_mass_vector,
           weight_mass_3tensor
    
    export IntegrationMethod, ExactIntegration, QuadratureIntegration

    abstract type AbstractGenerator{T} end
    @inline Base.eltype(::AbstractGenerator{T}) where T = T
    @inline Base.length(gen::AbstractGenerator) = gen.size
    @inline getpolynomial(gen::AbstractGenerator, n::Int) = gen[n]
    @inline getderivpolynomial(gen::AbstractGenerator, n::Int) = getderivpolynomial(gen)[n]

    include("fem/generators.jl")
    include("fem/basis.jl")
    include("fem/computations.jl")
    include("fem/utils.jl")
    include("fem/matrices.jl")
    include("fem/integration_formula.jl")


    # KOHN-SHAM MODEL
    export ExchangeCorrelation,NoExchangeCorrelation, SlaterXα, LSDA
    export exc, vxc, vxcUP, vxcDOWN
    export KohnShamExtended, ReducedHartreeFock
    include("models.jl")
    
    # SOLVER &CO
    abstract type KohnShamDiscretization end

    abstract type SCFMethod end
    abstract type SCFCache end
    abstract type SCFSolution end



    export LogConfig, LogBook
    include("log.jl")

    export KohnShamSolver, SolverOptions
    include("solver.jl")


    loopheader!(::SCFCache, ::SCFMethod, ::KohnShamSolver)      = nothing
    performstep!(::SCFCache, m::SCFMethod, ::KohnShamSolver)    = @warn "No performstep for the method $(typeof(m))"
    loopfooter!(::SCFCache, ::SCFMethod, ::KohnShamSolver)      = nothing
    monitor(::SCFCache, ::SCFMethod, ::KohnShamSolver)          = nothing
    register!(::SCFCache, ::SCFMethod, ::KohnShamSolver)        = nothing
    create_cache_method(m::SCFMethod, 
                        ::KohnShamDiscretization)               = @warn "No creation of cache for the method $(typeof(m))"
    switch!(c2::SCFCache, c1::SCFCache)                         = @warn "No way to init $(typeof(c2)) from of cache for the method ($(typeof(c1))"
    
    # DISCRETIZATION
    export LDADiscretization, LSDADiscretization
    include("discretization/lda.jl")
    include("discretization/lsda.jl")

    export DFTProblem
    include("problem.jl")

    ## SCF METHODS
    export CDA, ODA, Quadratic
    include("methods/rca.jl")
    include("methods/oda_procedure.jl")
    #include("methods/quadratic.jl")
    include("methods/combined.jl")

    export aufbau!
    include("aufbau.jl")

    export KohnShamSolution, eigenvector, density, total_charge
    include("solution.jl")

    export groundstate
    include("groundstate.jl")
end
