# STRUCTURE ATOM PROBLEM
struct AtomProblem
    T           # Type of numbers
    lh          # Truncation of orbital for l
    nh          # Truncation of orbital for n
    method      # SCF Method
    model       # Model
    Rmax        # Spatial cut-off
    Nmesh       # Number of points of the discretization
    typemesh    # Type of the Mesh
    optsmesh    # Options for the Mesh
    typebasis   # Type of the basis
    optsbasis   # Options for the basis
    typediscre  # Type of discretization
    name        # Name of the problem
    solveropts  # Option for the solvers

    AtomProblem(;T, lh, nh, method, model, Rmax, Nmesh,typemesh, optsmesh, typebasis, optsbasis, typediscre, name = "", kwargs...) = 
        new(T, lh, nh, method, model, Rmax, Nmesh,typemesh, optsmesh, typebasis, optsbasis, typediscre, name, kwargs)

    function AtomProblem(prob; 
        T = prob.T, lh = prob.lh, nh = prob.nh, method = prob.method, model = prob.model, Rmax = prob.Rmax, Nmesh = prob.Nmesh,
        typemesh = prob.typemesh, optsmesh = prob.optsmesh, typebasis = prob.typebasis, 
        optsbasis = prob.optsbasis, typediscre = prob.typediscre, name = prob.name, kwargs = prob.solveropts) 
        new(T, lh, nh, method, model, Rmax, Nmesh,typemesh, optsmesh, typebasis, optsbasis, typediscre, name, kwargs)
    end
end

function KohnShamResolution.groundstate(prob::AtomProblem)
    mesh = prob.typemesh(zero(prob.T), prob.Rmax, prob.Nmesh; T = prob.T, prob.optsmesh...)
    basis = prob.typebasis(mesh, prob.T; prob.optsbasis...)
    discretization = prob.typediscre(prob.lh, basis, mesh, prob.nh)
    groundstate(prob.model, discretization, prob.method; name = prob.name, prob.solveropts...)
end


# STRUCTURE ATOM FOR ANALYSING CONVERGENCE
struct AtomConvergenceNmesh
    vecNmesh        # Set of Nmesh used
    Error           # Dict of errors on orbitals energy : for each problem,
                    # there is a vector of errors depending on Nmesh
    num             # Number of orbitals energy used 
end

struct AtomConvergenceRmax
    vecRmax         # Set of Nmesh used
    Error           # Dict of errors on orbitals energy : for each problem,
                    # there is a vector of errors depending on Rmax
    num             # Number of orbitals energy used 
end