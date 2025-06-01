#--------------------------------------------------------------------
#                           SOLVER OPTIONS
#--------------------------------------------------------------------
struct SolverOptions{T <: Real, 
                    intexcType <: IntegrationMethod, 
                    intfemType <: IntegrationMethod}
    scftol::T                               # SCF tolerance
    maxiter::Int                            # Maximum of iteration done

    exc_integration_method::intexcType      # Integration method to compute integrals 
                                            # involving exchange correlation
    fem_integration_method::intfemType      # Integration method to compute integrals
                                            # for fem's matrices

    degen_tol::T                            # Tolerance to consider degenescence of orbitals energy

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
                exc_integration_method::IntegrationMethod = ExactIntegration(),
                fem_integration_method::IntegrationMethod = ExactIntegration(),
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
- `exc_integration_method`: Integration method for exchange–correlation terms.
- `fem_integration_method`: Integration method for FEM matrix assembly.
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
mutable struct KSESolver{   discretizationType <: KSEDiscretization,
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
    energies::Dict{Symbol,dataType}             # Storages of the energies
    logbook::logbookType                        # LogBook
    
    function KohnShamSolver(model::KSEModel, 
                            discretization::KSEDiscretization, 
                            alg::SCFAlgorithm; 
                            scftol::Real, 
                            maxiter::Int = 100,
                            exc_integration_method::IntegrationMethod = ExactIntegration(),
                            fem_integration_method::IntegrationMethod = ExactIntegration(),
                            degen_tol::Real = eps(eltype(discretization.basis)),
                            logconfig = LogConfig(),
                            verbose::Int = 3)
    
        # Set the data type as the one of the discretization basis
        T = discretization.elT
        
        # Init Cache of the Discretisation
        init_cache!(discretization, model, hartree, fem_integration_method)
        
        # Init Cache of the Method
        cache = create_cache_alg(alg, discretization)
        
        # Init Energies 
        energies = zero_energies(discretization, model)
        
        #  SolverOptions
        opts = SolverOptions(   T(scftol), 
                                maxiter, 
                                exc_integration_method, 
                                fem_integration_method, 
                                T(degen_tol),
                                UInt8(verbose))
        
        # Init log parameters
        niter = 0
        stopping_criteria = zero(T)
        
        logbook = LogBook(logconfig, T)

        new{typeof(discretization),
            typeof(model),
            typeof(alg),
            typeof(cache),
            typeof(opts),
            typeof(logbook),
            typeof(stopping_criteria)}( niter, 
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