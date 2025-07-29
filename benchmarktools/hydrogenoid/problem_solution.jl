# STRUCTURE HYDROGENOID PROBLEM

struct HydrogenoidOption
    nλ::Any         # Number of the eigenvalue wanted
    nU::Any         # Number of the eigenvector wanted
end

struct HydrogenoidProblem
    T::Any           # Type of numbers
    z::Any           # Number of electrons
    l::Any           # 2nd quantum number of orbital
    Rmax::Any        # Spatial cut-off
    Nmesh::Any       # Number of points of the discretization
    typemesh::Any    # Type of the Mesh
    optsmesh::Any    # Options for the Mesh
    typebasis::Any   # Type of the basis
    optsbasis::Any   # Options for the basis
    opts::Any        # Option for for the solution
    name::Any        # Name of the problem

    function HydrogenoidProblem(;
            T, z, l, Rmax, Nmesh, typemesh, typebasis, optsmesh, optsbasis,
            nλ = 1:10, nU = nλ, name = "")
        _nU = isnothing(nU) ? (1:0) : nU
        opts = HydrogenoidOption(nλ, _nU)
        new(T, z, l, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis, opts, name)
    end

    function HydrogenoidProblem(prob;
            T = prob.T, z = prob.z, l = prob.l, Rmax = prob.Rmax, Nmesh = prob.Nmesh,
            typemesh = prob.typemesh, typebasis = prob.typebasis,
            optsmesh = prob.optsmesh, optsbasis = prob.optsbasis,
            nλ = prob.nλ, nU = prob.nU, name = prob.name)
        _nU = isnothing(nU) ? (1:0) : nU
        opts = HydrogenoidOption(nλ, _nU)
        new(T, z, l, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis, opts, name)
    end
end

# STRUCTURE HYDROGENOID SOLUTION

struct HydrogenoidSolution
    problem::Any         # Original problem
    nλ::Any              # Number of the eigenvalue computed
    nU::Any              # Number of the eigenvector computed
    λ::Any               # Numerical eigenvalue
    U::Any               # Numerical eigenvector
    λtheo::Any           # Theoretical eigenvalue
    Utheo::Any           # Theoretical eigenvector
    Δλ::Any              # |λ - λtheo|
    ΔU::Any              # norm(U - Utheo) TO CLARIFY THE COMPUTATION OF THE NORM

    function HydrogenoidSolution(problem, λ, U)
        nλ = problem.opts.nλ
        nU = problem.opts.nU
        _λ = λ[nλ]
        _U = nothing
        λtheo = theoretical_eigenvalue(problem)
        Utheo = theoretical_eigenvector(problem)
        Δλ = abs.(_λ .- λtheo)
        ΔU = Δλ      ## TO DO
        new(problem, nλ, nU, _λ, _U, λtheo, Utheo, Δλ, ΔU)
    end
end

# STRUCTURE HYDROGENOID FOR ANALYSING CONVERGENCE

struct HydrogenoidConvergenceNmesh
    probs::Any           # Set of problems
    vecNmesh::Any        # Set of Nmesh used
    ΔΛ::Any              # Dict of errors on eigenvalues : for each problem,
    # there is a vector of errors depending on Nmesh
    ΔU::Any              # Dict of errors on eigenvectors : for each problem,
    # there is a vector of errors depending on Nmesh
    num::Any             # Number of eigenvalue and eigenvector used 
end

struct HydrogenoidConvergenceRmax
    probs::Any           # Set of problems
    vecRmax::Any         # Set of Rmax used
    ΔΛ::Any              # Dict of errors on eigenvalues : for each problem,
    # there is a vector of errors depending on Rmax
    ΔU::Any              # Dict of errors on eigenvectors : for each problem,
    # there is a vector of errors depending on Rmax
    num::Any             # Number of eigenvalue and eigenvector used
end
