
#  MAIN FUNCTION : GROUNDSTATE 
function groundstate(model::AbstractDFTModel, 
                     discretization::KohnShamDiscretization, 
                     method::SCFMethod; 
                     name::String = "", kwargs...)
    solver = init(model, discretization, method; kwargs...)
    solve!(solver)
    makesolution(solver, name)
end

groundstate(problem::DFTProblem; kwargs...) =  groundstate(model(problem), discretization(problem), method(problem); kwargs...)

#  SOLVE
function solve!(solver::KohnShamSolver)
    while (solver.stopping_criteria > solver.opts.scftol || iszero(solver.niter)) && solver.niter < solver.opts.maxiter
        loopheader!(solver)
        performstep!(solver.cache, solver.method, solver)
        loopfooter!(solver)
        monitor(solver)
    end
end

# INITALIZATION
function init(  model::AbstractDFTModel, discretization::KohnShamDiscretization, method::SCFMethod; 
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
    
    KohnShamSolver(niter, stopping_criteria, discretization, model, method, cache, opts, energies, logbook)
end


# LOOPHEADER
function loopheader!(solver::KohnShamSolver)
    # LOOPHEADER SPECIFIC FOR THE METHOD 
    loopheader!(solver.cache, solver.method, solver)
    nothing                                         
end 


# LOOPFOOTER
function loopfooter!(solver::KohnShamSolver)
    # INCREASE THE NUMBER OF ITERATIONS DONE
    solver.niter += 1

    # COMPUTE THE NEW STOPPING CRITERIA   
    stopping_criteria!(solver) 

    # UPDATE THE LOG                                           
    register!(solver)     
    
    # LOOPFOOTER SPECIFIC FOR THE METHOD
    loopfooter!(solver.cache, solver.method, solver) 
    nothing                                            
end 


# MONITOR : DISPLAY CURRENT STATE OF SOLVER
function monitor(solver::KohnShamSolver)
    if solver.opts.verbose > 0
        println("--------------------------")
        println("Iteration : $(solver.niter)")
    end
    if solver.opts.verbose > 1
        println("Selected Method : $(name(solver.method))")
        println("Stopping criteria: $(solver.stopping_criteria)")
        println("Total Energy: $(solver.energies[:Etot])")
    end
    if solver.opts.verbose > 2
        monitor(solver.cache, solver.method, solver)
    end
end


# COMPUTE THE STOPPING CRITERIA
function stopping_criteria!(solver::KohnShamSolver)
    solver.stopping_criteria =  stopping_criteria!(solver.cache, solver.method, solver)
end
    


# UPDATE THE LOG : STORE INTERMEDIATE STATE
function register!(solver::KohnShamSolver)
    @unpack logbook = solver
    @unpack occupation_number, orbitals_energy, stopping_criteria, energy, density = logbook.config

    # STORE THE STOPPING CRITERIA 
    !stopping_criteria || push!(solver.logbook.stopping_criteria_log, solver.stopping_criteria)   
    
    #=
    # STORE THE ORBITALS ENERGY
    orbitals_energy ?   push!(solver.logbook.orbitals_energy_log, copy(solver.ϵ))              : nothing     
    
    # STORE THE OCCUPATION NUMBERS
    occupation_number ? push!(solver.logbook.occupation_number_log, copy(solver.n))          : nothing    
    =#

    # STORE THE TOTAL ENERGY
    !energy || push!(solver.logbook.energy_log, solver.energies[:Etot])

    # REGISTER DATA SPECIFIC TO THE METHODS
    register!(solver.cache, solver.method, solver)
end



# MAKE THE SOLUTION
function makesolution(solver::KohnShamSolver, name::String)
    datas = makesolution(solver.cache, solver.method, solver)
    KohnShamSolution(solver, name, datas)
end