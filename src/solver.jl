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
    
    function KohnShamSolver(niter::Int, 
                            stopping_criteria::Real, 
                            discretization::KohnShamDiscretization, 
                            model::AbstractDFTModel, 
                            method::SCFMethod,
                            cache::SCFCache,
                            opts::SolverOptions,  
                            energies::Dict{Symbol,<:Real}, 
                            logbook::LogBook)
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