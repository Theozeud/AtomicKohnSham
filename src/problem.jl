#--------------------------------------------------------------------
#                         PROBLEM STRUCTURE
#--------------------------------------------------------------------
"""
    AtomProblem(; T, lh, nh, alg, model, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis,
                  typediscre, name = "", kwargs...)

Defines an atomic simulation problem.

This structure gathers all the information required to set up and solve the electronic ground state of an atom or ion using a 
self-consistent field (SCF) procedure, with control over the discretization, mesh, basis, and solver options.

# Keyword Arguments
- `T`: Number type used in computations (e.g., `Float64` or `DoubleFloats`).
- `lh::Int`: Maximum value of the orbital angular momentum quantum number `l`.
- `nh::Int`: Maximum principal quantum number `n`.
- `alg`: SCF method (e.g., `ODA`, `Quadratic`).
- `model`: The physical model to be solved (`KSEModel`).
- `Rmax::T`: Radial cut-off, i.e., the domain is `[0, Rmax]`.
- `Nmesh::Int`: Number of mesh points for discretization.
- `typemesh`: Type of mesh used (e.g., `ExponentialMesh`).
- `optsmesh`: Options passed to the mesh constructor.
- `typebasis`: Type of finite element basis (e.g., `IntLegendre`).
- `optsbasis`: Options passed to the basis constructor.
- `name::String`: Optional name of the problem (default: empty string).
- `kwargs`: Additional keyword arguments forwarded to the solver configuration.
"""
struct AtomProblem{ T <: Real,
                    A <: SCFAlgorithm,
                    M <: KSEModel,
                    S <: AbstractString,
                    OM <: NamedTuple,
                    OB <: NamedTuple,
                    meshType <: Mesh,
                    basisType <: FEMBasis
                    }
    lh::Int             # Truncation of orbital for l
    nh::Int             # Truncation of orbital for n
    alg::A              # SCF Algorithm
    model::M            # Model
    Rmax::T             # Spatial cut-off
    Nmesh::Int          # Number of points of the discretization
    optsmesh::OM        # Options for the Mesh
    optsbasis::OB       # Options for the basis
    name::S             # Name of the problem
    solveropts::OS      # Option for the solvers

    function AtomProblem(;T, lh, nh, alg, model, Rmax, Nmesh,typemesh, optsmesh, typebasis, optsbasis, name = "", kwargs...)
        new{T, typeof(alg), typeof(model), typeof(name), typeof(optsmesh),
            typeof(optsbasis), typemesh, typebasis}(lh, nh, alg, model, Rmax, Nmesh, optsmesh, optsbasis, name, kwargs)
    end

    function AtomProblem(prob; 
        T = prob.T, lh = prob.lh, nh = prob.nh, alg = prob.alg, model = prob.model, Rmax = prob.Rmax, Nmesh = prob.Nmesh,
        typemesh = prob.typemesh, optsmesh = prob.optsmesh, typebasis = prob.typebasis, 
        optsbasis = prob.optsbasis, name = prob.name, kwargs = prob.solveropts) 
        new{T, typeof(alg), typeof(model), typeof(name), typeof(optsmesh),
            typeof(optsbasis), typemesh, typebasis}(lh, nh, alg, model, Rmax, Nmesh, optsmesh, optsbasis, name, kwargs)
    end
end