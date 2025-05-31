using KohnShamResolution
using TimerOutputs
using UnPack
using TensorOperations
using LinearAlgebra

using KohnShamResolution:   init, loopheader!, loopfooter!, makesolution,
                            prepare_eigenvalue_problem!, find_orbital!, aufbau!, density!, update_density!,
                            compute_total_energy, compute_kinetic_energy, compute_coulomb_energy, compute_hartree_energy, isthereExchangeCorrelation, compute_exchangecorrelation_energy,
                            tensor_matrix_dict!, compute_density,
                            init_cache!, create_cache_method, init_energies, SolverOptions



to = TimerOutput()


#####################################################################
#                               PARAMETERS
#####################################################################
model = ReducedHartreeFock(1, 1)
method = ODA(0.3)
lh = 0
nh = 1
Nh = 20
Rmin = 0
Rmax = 100
ordermax = 20

#####################################################################
#                           INITIALIZATION
#####################################################################

@timeit to "Create mesh" m = expmesh(Rmin, Rmax, Nh; s = 1.5)

@timeit to "Create basis" basis = P1IntLegendreGenerator(m; ordermax = ordermax)

@timeit to "init Discretization" discretization = LDADiscretization(lh, basis, m, nh)

@timeit to "Init Solver" begin 

        # Set the data type as the one of the discretization basis
        T = discretization.elT

        # Init Cache of the Discretisation
        @timeit to "Init Cache discretization" init_cache!(discretization, model, true, ExactIntegration())
    
        # Init Cache of the Method
        @timeit to "Create Cache Method" cache = create_cache_method(method, discretization)
    
        # Init Energies 
        @timeit to "Init Energies" energies = init_energies(discretization, model)
          
        #  SolverOptions
        @timeit to "Solver Options" opts = SolverOptions(T(1e-10), 
                                                        100, 
                                                        ExactIntegration(), 
                                                        ExactIntegration(), 
                                                        T(true), 
                                                        eps(T),
                                                        UInt8(0))
    
        # Init log parameters
        niter = 0
        stopping_criteria = zero(T)
        
        @timeit to "Create logbook" logbook = LogBook(LogConfig(), T)
        
        solver = KohnShamSolver(niter, stopping_criteria, discretization, model, method, cache, opts, energies, logbook)
end


#####################################################################
#                               SOLVE
#####################################################################

for i ∈ 1:10

    @timeit to "Loop header" loopheader!(solver)


    @timeit to "PerformStep" begin

        @unpack D, Dprev, U, ϵ, n = cache
    
        # STEP 1 : PREPARE THE EIGENVALUE PROBLEM
        @timeit to "Prepare eigenvalue problem" prepare_eigenvalue_problem!(discretization, model, Dprev, opts.hartree)

        # STEP 2 : FIND ORBITALS AND CORRESPONFING ENERGIES
        @timeit to "Find orbital" find_orbital!(discretization, U, ϵ)

        # STEP 3 : FILL THE OCCUPATION NUMBER MATRIX ACCORDINGLY WITH THE AUFBAU PRINCIPLE
        @timeit to "aufbau" aufbau!(cache, solver)

        if !cache.flag_degen

            # STEP 4 : COMPUTE A GUESS DENSITY
            @timeit to "density computation" density!(discretization, U, n, D)

            # STEP 5 : COMPUTE ALL ENERGIES
            @timeit to "compute energy" begin 
                @timeit to "Etot" energies[:Etot] = compute_total_energy(discretization, model, D, n, ϵ)
                @timeit to "Ekin" energies[:Ekin] = compute_kinetic_energy(discretization, U, n)
                @timeit to "Ecou" energies[:Ecou] = compute_coulomb_energy(discretization, U, n)
                @timeit to "Ehar" energies[:Ehar] = compute_hartree_energy(discretization, D)
            end
        end

      # STEP  6 : COMPUTE THE NEW DENSITY
      @timeit to "update density" update_density!(cache, method, solver)  
    end

    @timeit to "Loop footer" loopfooter!(solver)
end



@timeit to "Make Solution" makesolution(solver, "")




#####################################################################
#                    PRINT TIME DATAS IN A FILE
#####################################################################

original_stdout = stdout
output_file = open("tests/Performance/solver/oda-rHF-hydro.txt", "a")
redirect_stdout(output_file)

println()
println("==============================================")
println("PERFORMANCE PARAMETERS")
println("Nh = $Nh")
println("lh = $lh")
println("ordermax = $ordermax")
println("Method = ODA")
println("Number of basis function : $(length(basis))")
println("==============================================")

print(to)

redirect_stdout(original_stdout)
close(output_file)
println("Performance finished")
