struct SolverOptions{T <: Real, 
                    intexcType <: IntegrationMethod, 
                    intfemType <: IntegrationMethod}
    scftol::T                               # SCF tolerance
    maxiter::Int                            # Maximum of iteration done

    exc_integration_method::intexcType      # Integration method to compute integrals 
                                            # involving exchange correlation
    fem_integration_method::intfemType      # Integration method to compute integrals
                                            # for fem's matrices

    hartree::T                              # Coefficient multiply to the Hartree Matrix : 
                                            # 0 -> no hartree term, 
                                            # 1-> full hartree term

    degen_tol::T                            # Tolerance to consider degenescence of orbitals energy

    verbose::UInt8                          # Say how many details are printed at the end 
                                            # of each iterations :
                                            # 0 : Zero details
                                            # 1 : Iterations
                                            # 2 : Computations                
end


mutable struct KohnShamSolver{  discretizationType <: KohnShamDiscretization,
                                modelType <: AbstractDFTModel,
                                methodType <: SCFMethod,
                                cacheType <: SCFCache,
                                optsType <: SolverOptions,
                                logbookType <: LogBook,
                                dataType <: Real}
                            
    niter::Int                                  # Number of iterations
    stopping_criteria::dataType                 # Current stopping criteria
    discretization::discretizationType          # Discretization parameters
    model::modelType                            # Model
    method::methodType                          # Iterative method
    cache::cacheType                            # Cache associated to the method
    opts::optsType                              # Solver options
    energies::Dict{Symbol,dataType}             # Storages of the energies
    logbook::logbookType                        # LogBook
    
    function KohnShamSolver(model::AbstractDFTModel, 
                            discretization::KohnShamDiscretization, 
                            method::SCFMethod; 
                            scftol::Real, 
                            maxiter::Int = 100,
                            exc_integration_method::IntegrationMethod = ExactIntegration(),
                            fem_integration_method::IntegrationMethod = ExactIntegration(),
                            hartree::Real = 1, 
                            degen_tol::Real = eps(eltype(discretization.basis)),
                            logconfig = LogConfig(),
                            verbose::Int = 3)
    
        # Set the data type as the one of the discretization basis
        T = discretization.elT
        
        # Init Cache of the Discretisation
        init_cache!(discretization, model, hartree, fem_integration_method)
        
        # Init Cache of the Method
        cache = create_cache_method(method, discretization)
        
        # Init Energies 
        energies = init_energies(discretization, model)
        
        #  SolverOptions
        opts = SolverOptions(T(scftol), 
                    maxiter, 
                    exc_integration_method, 
                    fem_integration_method, 
                    T(hartree), 
                    T(degen_tol),
                    UInt8(verbose))
        
        # Init log parameters
        niter = 0
        stopping_criteria = zero(T)
        
        logbook = LogBook(logconfig, T)

        new{typeof(discretization),
            typeof(model),
            typeof(method),
            typeof(cache),
            typeof(opts),
            typeof(logbook),
            typeof(stopping_criteria)}( niter, 
                                        stopping_criteria, 
                                        discretization, 
                                        model, 
                                        method, 
                                        cache, 
                                        opts,
                                        energies,
                                        logbook)
    end
end