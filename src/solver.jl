#--------------------------------------------------------------------
#                           SOLVER OPTIONS
#--------------------------------------------------------------------
struct SolverOptions{T <: Real}
    scftol::T                               # SCF tolerance
    maxiter::Int                            # Maximum of iteration done

    degen_tol::T                            # Tolerance to consider degenescence of orbitals
                                            # energy

    verbose::UInt8                          # Say how many details are printed at the end
    # of each iterations :
    # 0 : Zero details
    # 1 : Iterations
    # 2 : Computations
end

#--------------------------------------------------------------------
#                           SOLVER STRUCTURE
#--------------------------------------------------------------------
"""
    KSESolver(  model::KSEModel,
                discretization::KSEDiscretization,
                alg::SCFAlgorithm;
                scftol::Real,
                maxiter::Int = 100,
                degen_tol::Real = eps(eltype(discretization.basis)),
                logconfig = LogConfig(),
                verbose::Int = 3) -> KSESolver

Builds a Kohn-Sham-Extended solver object to apply the self-consistent field (SCF) algorithms.

The `KSESolver` contains all internal data required for iterative SCF procedures including the model, discretization, algorithm, and solver options.
Use `solve!(solver)` to run the self-consistent procedure, and access solution data via `KSESolution(solver)`.

# Arguments
- `model`: The Kohn–Sham extended model (e.g. `RHF`, `SlaterXα`).
- `discretization`: Spatial discretization structure.
- `alg`: SCF algorithm to use (e.g. `ODA``, `Quadratic()`).
- `scftol`: Tolerance on the SCF convergence.
- `maxiter`: Maximum number of SCF iterations.
- `degen_tol`: Tolerance to detect degeneracy in eigenvalues (default: machine epsilon).
- `logconfig`: Logging configuration (`LogConfig()`).
- `verbose`: Verbosity level (0: silent, 3: detailed).

# Fields
- `niter::Int`: Current number of iterations.
- `stopping_criteria::Real`: Value of convergence metric (e.g., residual density norm).
- `discretization::KSEDiscretization`: Discretization structure (mesh, basis, etc.).
- `model::KSEModel`: The model (e.g., RHF, Slater).
- `alg::SCFAlgorithm`: The SCF algorithm used.
- `cache::SCFCache`: Internal cache for accelerating SCF iteration.
- `opts::SolverOptions`: Options related to tolerances, methods, and iteration limits.
- `energies::Dict{Symbol,Real}`: Dictionary storing energy contributions (total, kinetic, Hartree...).
- `logbook::LogBook`: Structure controlling what is logged and storing convergence history.

# Returns
- A fully initialized `KSESolver` ready for use with `solve!`.
"""
mutable struct KSESolver{discretizationType <: KSEDiscretization,
    modelType <: KSEModel,
    algType <: SCFAlgorithm,
    cacheType <: SCFCache,
    optsType <: SolverOptions,
    logbookType <: LogBook,
    dataType <: Real}
    niter::Int                                  # Number of iterations
    stopping_criteria::dataType                 # Current stopping criteria
    discretization::discretizationType          # Discretization parameters
    model::modelType                            # Model
    alg::algType                                # SCF Algorithm
    cache::cacheType                            # Cache associated to the algorithm
    opts::optsType                              # Solver options
    energies::Dict{Symbol, dataType}             # Storages of the energies
    logbook::logbookType                        # LogBook

    function KSESolver(model::KSEModel,
            discretization::KSEDiscretization,
            alg::SCFAlgorithm;
            scftol::Real,
            maxiter::Int = 100,
            degen_tol::Real = eps(eltype(discretization.basis)),
            logconfig = LogConfig(),
            verbose::Int = 0)

        # Set the data type as the one of the discretization basis
        T = eltype(discretization)

        # Init Cache of the Discretisation
        init_cache!(discretization, model)

        # Init Cache of the Method
        cache = create_cache_alg(alg, discretization)

        # Init Energies
        energies = zero_energies(discretization, model)

        #  SolverOptions
        opts = SolverOptions(T(scftol),
            maxiter,
            T(degen_tol),
            UInt8(verbose))

        # Init log parameters
        niter = 0
        stopping_criteria = zero(T)

        logbook = LogBook(logconfig, T, energies, alg)

        new{typeof(discretization),
            typeof(model),
            typeof(alg),
            typeof(cache),
            typeof(opts),
            typeof(logbook),
            typeof(stopping_criteria)}(niter,
            stopping_criteria,
            discretization,
            model,
            alg,
            cache,
            opts,
            energies,
            logbook)
    end
end
