function scf_step!(alg::ODA, solver::KSESolver)
    @unpack discretization, model, energies, algcache, niter = solver
    @unpack D, Dprev, U, ϵ, n, aufbaucache, energies_prev  = algcache

    # ==== STEP 1 : Save previous density and energies
    @. Dprev = D
    energies_prev.Etot = energies.Etot
    energies_prev.Ekin = energies.Ekin
    energies_prev.Ecou = energies.Ecou
    energies_prev.Ehar = energies.Ehar
    energies_prev.Eexc = energies.Eexc

    # ==== STEP 2 : Build the hamiltonian
    assemble_hamiltonian!(discretization, model, Dprev)

    # ==== STEP 3 : Find the orbitals and associated energy
    find_orbital!(discretization, U, ϵ)

    # ==== STEP 4 : Fill the occupation numbers accordingly the aufbau method
    aufbau!(n, ϵ, U, model, discretization, aufbaucache, niter)

    # ==== STEP 5 : Compute the trial density
    if aufbaucache.postcomputations # VARIABLE SHOULD BE IN ODA CACHE
        D .= aufbaucache.D1
        energies.Etot = aufbaucache.energies.Etot
        energies.Ekin = aufbaucache.energies.Ekin
        energies.Ecou = aufbaucache.energies.Ecou
        energies.Ehar = aufbaucache.energies.Ehar
        energies.Eexc = aufbaucache.energies.Eexc
    else
        density!(discretization, U, n, D)
        energies.Etot = compute_total_energy(discretization, model, D, n, ϵ)
        energies.Ekin = compute_kinetic_energy(discretization, U, n)
        energies.Ecou = compute_coulomb_energy(discretization, U, n)
        energies.Ehar = compute_hartree_energy(discretization, D)
        energies.Eexc = compute_exc_energy(discretization, model, D)
    end

    # ==== STEP 6 : Relaxation Step
    niter > 0 ? line_search_density!(alg, solver) : nothing

    # ==== STEP 7 : Return the current stopping criteria
    stop_D = norm(D - Dprev)
    stop_Etot = abs(energies.Etot - energies_prev.Etot)
    return max(stop_D, stop_Etot)
end


scf_converged(alg::ODA, solver::KSESolver) = solver.stopping_criteria < alg.scftol



function monitor(cache::ODACache, ::ODA, ::KSESolver)
    println("Relaxed Parameter : $(cache.t)")
end

function postcomputations!(algcache::ODACache, solver::KSESolver)
    @unpack discretization, model = solver
    assemble_hartree_pot!(discretization, algcache.D; coeff = model.hartree)
end

"""
    find_orbital!(disc, U, ϵ)

Solve the discretized Kohn–Sham eigenvalue problems.

For each angular momentum channel and spin component, this function solves
the generalized eigenvalue problem and stores the lowest eigenvalues in `ϵ`
and the corresponding orbitals in `U`.
"""
function find_orbital!(discretization::KSEDiscretization,
                      U::AbstractArray{<:Real}, ϵ::AbstractArray{<:Real})
    @unpack lₕ, nₕ, nspin, ksham, femops = discretization
    # Solve the generalized eigenvalue problem for each section (l,σ)
    for σ in 1:nspin
        for l in 0:lₕ
            @views vH = ksham.H[:, :, l + 1, σ]
            λ, V = eigen(Symmetric(femops.S*vH*femops.S))
            @views ϵ[l + 1, :, σ] = λ[1:nₕ]
            @views U[:, :, l + 1, σ] = femops.S*V[:, 1:nₕ]
        end
    end
end

function line_search_density!(alg::ODA, solver::KSESolver)
    @unpack discretization, energies, algcache = solver
    @unpack D, Dprev, Dbuf, energies_prev, F = algcache
    @unpack frozen_t, maxiter_ls, abstol_ls, reltol_ls = alg

    energy_har01 = compute_hartree_mix_energy(discretization, D, Dprev)

    if frozen_t
        tnew = algcache.t
        # Update next density
        @. D = tnew * D + (1 - tnew) * Dprev
        # Compute the new energies
        energies.Ekin = tnew*energies.Ekin + (1-tnew)*energies_prev.Ekin
        energies.Ecou = tnew*energies.Ecou + (1-tnew)*energies_prev.Ecou
        energies.Ehar = (tnew^2)*energies.Ehar + ((1-tnew)^2)*energies_prev.Ehar +
                            2 * tnew * (1-tnew) * energy_har01
        energies.Eexc = F isa Real ? 0 : F(tnew)
        energies.Etot = energies.Ekin + energies.Ecou + energies.Ehar + energies.Eexc

    else
        # Find the optimal t∈[0,1]
        algcache.t, energies.Etot = line_search_energy(energies.Ekin, energies_prev.Ekin,
                                                   energies.Ecou, energies_prev.Ecou,
                                                   energies.Ehar, energies_prev.Ehar,
                                                   energy_har01, F; maxiter = maxiter_ls,
                                                   abstol = abstol_ls, reltol = reltol_ls)
        tnew = algcache.t
        # Update next density
        @. D = tnew * D + (1 - tnew) * Dprev
        # Compute the new energies
        energies.Ekin = tnew*energies.Ekin + (1-tnew)*energies_prev.Ekin
        energies.Ecou = tnew*energies.Ecou + (1-tnew)*energies_prev.Ecou
        energies.Ehar = tnew^2*energies.Ehar + (1-tnew)^2*energies_prev.Ehar +
                            2 * tnew * (1-tnew) * energy_har01
        energies.Eexc = energies.Etot - energies.Ekin - energies.Ecou - energies.Ehar
    end
    return nothing
end
