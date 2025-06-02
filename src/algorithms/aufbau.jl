# RETURN THE MAXIMAL NUMBER OF ELECTRONS THAT CAN BE PUT AT THE LAYER CORRESPONDING TO THE INDEX idx

function convert_index(discretization::KSEDiscretization, idx::Int)
    @unpack lₕ  = discretization
    l = rem(idx - 1, lₕ+1)
    k = div(idx-1,lₕ+1)+1
    return (l,k)
end

#=
function convert_index(discretization::LSDADiscretization, idx::Int)
    @unpack lₕ, nₕ  = discretization
    σ       = div(idx-1, (lₕ+1)*nₕ)+1
    idxσ    = rem(idx-1, (lₕ+1)*nₕ)
    l = rem(idxσ, lₕ+1)
    k = div(idxσ,lₕ+1)+1
    return (l,k,σ)
end
=#
function degeneracy(discretization::KSEDiscretization, idx::Int)
    l,_ = convert_index(discretization::KSEDiscretization, idx)
    return 4 * l + 2
end

#=
function degeneracy(discretization::LSDADiscretization, idx::Int)
    l,_ = convert_index(discretization, idx)
    return 2 * l + 1
end
=#

#####################################################################
#                        AUFBAU PRINCIPLE
#####################################################################


function aufbau!(cache::RCACache, solver::KSESolver)

    @unpack model, discretization, opts, energies = solver
    @unpack U, ϵ, n, Noccup, D, tmpD, tmpD2, index_aufbau, energies_prev = cache

    # INIT OCCUPATION NUMBER
    fill!(n, zero(eltype(n))) 
    fill!(Noccup, zero(eltype(Noccup)))

    # COLLECT THE INDICES OF THE SORTED ORBITAL ENERGIES
    index_aufbau .= sortperm(vec(ϵ))
   
    # LOOP TO FILL THE ORBITALS
    remain = model.N
    idx = 1
    while remain > 0 && idx ≤ length(index_aufbau)

        # FIND ALL THE ORBITALS WITH THE SAME ENERGY
        indices_degen = [index_aufbau[idx]]  
        idx += 1
        while idx ≤ length(index_aufbau) && abs(ϵ[index_aufbau[idx]] - ϵ[first(indices_degen)]) < opts.degen_tol && length(indices_degen) <2
            push!(indices_degen, index_aufbau[idx])
            idx += 1
        end

        # COMPUTE DEGENERACY
        degen = zeros(Int, length(indices_degen))
        for i ∈ eachindex(indices_degen)
            degen[i] = degeneracy(discretization, indices_degen[i])
        end
        total_degen = sum(degen)

        # ELECTRON DISTRIBUTIONS
        if remain - total_degen ≥ 0
            # IN THIS CASE WE CAN FILL ALL THE LAYERS WITH THE SAME ORBITAL ENERGY

            for i in eachindex(indices_degen)
                n[indices_degen[i]] = degen[i]
                normalization!(discretization, U, convert_index(discretization,i)...) 
            end
            remain -= total_degen
            Noccup[1] += length(indices_degen)

        elseif length(indices_degen) == 1
            # IN THIS CASE, WE FILL THE ORBITAL WITH THE REMAIN ELECTRONS
            Noccup[2] += 1

            n[first(indices_degen)] = remain
            normalization!(discretization, U, convert_index(discretization,first(indices_degen))...)
            break

        elseif length(indices_degen) == 2
            # IN THIS CASE WE NEED TO FILL THE LAYERS WITH THE OPTIMAL REPARTITION
            cache.flag_degen = true
            Noccup[2] += 2

            # NORMALIZATION OF COEFFICIENTS OF ORBITAL
            for i ∈ indices_degen
                normalization!(discretization, U, convert_index(discretization,i)...)
            end

            # COMPUTE ENERGIES FOR ONE EXTREMA 
            degen_first = degen[1]
            if degen_first ≥ remain
                n[indices_degen[1]] = remain
                n[indices_degen[2]] = zero(remain)
                n1_0 = remain
                n2_0 = zero(remain)
            else
                n[indices_degen[1]] = degen_first
                n[indices_degen[2]] = remain - degen_first
                n1_0 = degen_first
                n2_0 = remain - degen_first
            end
            
            density!(discretization, U, n, D)
            energy_kin0 = compute_kinetic_energy(discretization, U, n)
            energy_cou0 = compute_coulomb_energy(discretization, U, n)
            energy_har0 = compute_hartree_energy(discretization, D)

            # COMPUTE ENERGIES FOR THE OTHER EXTREMA
            degen_second = degen[2]
            if degen_second ≥ remain
                n[indices_degen[1]] = zero(remain)
                n[indices_degen[2]] = remain
                n1_1 = zero(remain)
                n2_1 = remain
            else
                n[indices_degen[1]] = degen_first
                n[indices_degen[2]] = remain - degen_first
                n1_1 = remain - degen_second
                n2_1 = degen_second
            end
            
            density!(discretization, U, n, tmpD)
            energy_kin1 = compute_kinetic_energy(discretization, U, n)
            energy_cou1 = compute_coulomb_energy(discretization, U, n)
            energy_har1 = compute_hartree_energy(discretization, tmpD)

            energy_har01 = compute_hartree_mix_energy(discretization, D, tmpD)
            energy_har10 = compute_hartree_mix_energy(discretization, tmpD, D)

            # FIND THE OPTIMUM OCCUPATION
            cache.tdegen, energies[:Etot] = find_minima_oda(energy_kin0, energy_kin1, 
                                                            energy_cou0, energy_cou1, 
                                                            energy_har0, energy_har1, 
                                                            energy_har01, energy_har10,
                                                            D, tmpD, tmpD2, model, discretization)

            # UPDATE THE OCCUPATION NUMBERS
            n[indices_degen[1]] = cache.tdegen * n1_0 + (1-cache.tdegen) * n1_1
            n[indices_degen[2]] = cache.tdegen * n2_0 + (1-cache.tdegen) * n2_1

            # UPDATE THE DENSITY
            @. D = cache.tdegen * D + (1 - cache.tdegen) * tmpD

            # UPDATE THE ENERGIES
            t = cache.tdegen
            energies[:Ekin] = t*energy_kin0 + (1-t)*energy_kin1
            energies[:Ecou] = t*energy_cou0 + (1-t)*energy_cou1
            energies[:Ehar] = t^2*energy_har0 + (1-t)^2*energy_har1 + t*(1-t) * (energy_har01 + energy_har10)
            if isthereExchangeCorrelation(model)
                energies[:Eexc] = compute_exchangecorrelation_energy(discretization, model, D)
            end

            break
        else
            @error("This case of degeneracy is not coded.")
            break
        end

    end
    Noccup[3] = length(ϵ) - Noccup[1] - Noccup[2] 
    nothing
end


