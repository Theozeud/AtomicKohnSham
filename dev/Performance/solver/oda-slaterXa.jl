using KohnShamResolution
using TimerOutputs
using UnPack
using TensorOperations
using LinearAlgebra

using KohnShamResolution: init, loopheader!, loopfooter!, makesolution,
                          prepare_eigenvalue_problem!, find_orbital!, aufbau!, density!,
                          update_density!,
                          compute_total_energy, compute_kinetic_energy,
                          compute_coulomb_energy, compute_hartree_energy,
                          isthereExchangeCorrelation, compute_exchangecorrelation_energy,
                          tensor_matrix_dict!, compute_density

to = TimerOutput()

#####################################################################
#                               PARAMETERS
#####################################################################
model = SlaterXα(8, 8)
method = ODA(0.3)
lh = 1
nh = 2
Nh = 10
Rmin = 0
Rmax = 150
ordermax = 10

#####################################################################
#                           INITIALIZATION
#####################################################################

@timeit to "Create mesh" m = expmesh(Rmin, Rmax, Nh; s = 0.9)

@timeit to "Create basis" basis = P1IntLegendreGenerator(m; ordermax = ordermax)

@timeit to "init Discretization" discretization = LDADiscretization(lh, basis, m, nh)

@timeit to "Init Solver" solver = KohnShamResolution.init(
    model, discretization, method; scftol = 1e-3, hartree = false,
    logconfig = LogConfig(orbitals_energy = true))

#####################################################################
#                               SOLVE
#####################################################################

for i in 1:10
    @timeit to "Loop header" loopheader!(solver)

    @timeit to "PerformStep" begin
        @unpack opts, energies, cache = solver
        @unpack D, Dprev, U, ϵ, n = cache

        # STEP 1 : PREPARE THE EIGENVALUE PROBLEM
        @timeit to "Prepare eigenvalue problem" prepare_eigenvalue_problem!(
            discretization, model, Dprev, opts.hartree)

        # STEP 2 : FIND ORBITALS AND CORRESPONFING ENERGIES
        @timeit to "Find orbital" find_orbital!(discretization, U, ϵ)

        # STEP 3 : FILL THE OCCUPATION NUMBER MATRIX ACCORDINGLY WITH THE AUFBAU PRINCIPLE
        @timeit to "aufbau" aufbau!(cache, solver)

        if !cache.flag_degen

            # STEP 4 : COMPUTE A GUESS DENSITY
            @timeit to "density computation" density!(discretization, U, n, D)

            # STEP 5 : COMPUTE ALL ENERGIES
            @timeit to "compute energy" begin
                @timeit to "Etot" energies[:Etot] = compute_total_energy(
                    discretization, model, D, n, ϵ)
                @timeit to "Ekin" energies[:Ekin] = compute_kinetic_energy(
                    discretization, U, n)
                @timeit to "Ecou" energies[:Ecou] = compute_coulomb_energy(
                    discretization, U, n)
                @timeit to "Ehar" energies[:Ehar] = compute_hartree_energy(
                    discretization, D)
                @timeit to "Eexc" begin
                    ρ(x::Float64) = compute_density(discretization, D, x)
                    f(x::Float64, p::Float64) = exc(model.exc, ρ(x)) * x^2
                    ρ(1.0)
                    f(1.0, 1.0)
                    @timeit to "Compute density" ρ(1.0)
                    @timeit to "Eval objectiv function" f(1.0, 1.0)
                    @timeit to "Eexc" energies[:Eexc] = compute_exchangecorrelation_energy(
                        discretization, model, D)
                end
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
output_file = open("tests/Performance/solver/oda-SlaterXa.txt", "a")
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
