module KohnShamResolution

    using LinearAlgebra
    using SparseArrays
    using TimerOutputs
    using ProgressMeter
    using UnPack
    
    include("utils.jl")
    include("tools/computation_tools.jl")

    export LaurentPolynomial
    export Monomial, deg, degmax, degmin, haslog, ismonomial, iszero
    export integrate!, integrate, deriv!, deriv, scalar_product, elag!
    include("basis/laurentpolynomial.jl")

    export LaurentPolynomialBasis
    export mass_matrix, weight_mass_matrix
    include("basis/laurentpolynomialbasis.jl")

    export OneDMesh
    export mesh, find_index
    include("mesh.jl")

    export PiecewiseLaurentPolynomial
    export get_support
    include("basis/piecewiselaurentpolynomial.jl")

    export HatFunctionP1
    include("basis/hat_functions.jl")

    export KohnShamExtended
    export exchcorr, charge, nbelec, potential
    include("model/model.jl")

    export KohnShamDiscretization, KohnShamSphericalDiscretization
    include("discretization/abstract_discretization.jl")
    include("discretization/spherical_discretization.jl")

    include("methods/abstractmethods.jl")

    export DFTProblem
    include("problem.jl")

    include("solver.jl")

    export ODA
    include("methods/oda.jl")

    include("solution.jl")

    export groundstate
    include("groundstate.jl")
    

end
