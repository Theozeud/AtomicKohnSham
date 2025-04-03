#####################################################################
#                          RCA METHOD
#####################################################################

abstract type RCAMethod <: SCFMethod end

#####################################################################
#                          RCA CACHE
#####################################################################

mutable struct RCACache{dataType <: Real,
                        densityType <: AbstractArray,
                        orbitalsType <: AbstractArray,
                        orbitalsenergyType <: AbstractArray,
                        occupationType <: AbstractArray} <: SCFCache

    t::dataType                                 # Relaxation Parameter

    D::densityType                              # Density Matrix at current time
    Dprev::densityType                          # Density Matrix at previous time
    tmpD::densityType                           # Storage for extra Density Matrix 
                                                # (usefull for degeneracy)

    U::orbitalsType                             # Coefficient of orbitals at current time
    ϵ::orbitalsenergyType                       # Orbitals energy at current time
    n::occupationType                           # Occupation number at current time 
    Noccup::Vector{Int}                         # Triplet (Nf,Np,Nv) where  
                                                #     - Nf is the number of fully occupied  orbitals, 
                                                #     - Np is the number of partially occupied orbitals,
                                                #     - Nv is the number of virtual orbitals
    
    flag_degen::Bool                            # Flag to know if there is a degeneracy 
    tdegen::dataType                            # Interpolation parameters in case of degeneracy

    index_aufbau::Vector{Int}                   # Indices List for aufbau
    energies_prev::Dict{Symbol,dataType}        # Energies at previous time                                         
end

function create_cache_method(method::RCAMethod, discretization::KohnShamDiscretization)
    @unpack elT = discretization
    t                   = elT(method.t)
    D                   = init_density(discretization)
    Dprev               = init_density(discretization)
    tmpD                = init_density(discretization)
    U                   = init_orbitals(discretization)
    ϵ                   = init_orbitals_energy(discretization)
    n                   = init_occupation_number(discretization)
    Noccup              = zeros(Int,3)
    tdegen              = zero(elT)
    index_aufbau        = zeros(Int, length(ϵ))
    energies_prev       = Dict( :Etot => zero(elT),
                                :Ekin => zero(elT),
                                :Ecou => zero(elT),
                                :Ehar => zero(elT))
    RCACache{
        elT,
        typeof(D), 
        typeof(U), 
        typeof(ϵ), 
        typeof(n)}(t, D, Dprev, tmpD, U, ϵ, n, Noccup, false, tdegen, index_aufbau, energies_prev)
end


#####################################################################
#                          RCA SOLUTION
#####################################################################

struct RCASolution{ densityType <: AbstractArray,
                    orbitalsType <: AbstractArray,
                    orbitalsenergyType <: AbstractArray,
                    occupationType} <: SCFSolution
    density_coeffs::densityType                                 # Density Matrix at final time
    orbitals::orbitalsType                                      # Coefficient of orbitals at final time
    orbitals_energy::orbitalsenergyType                         # Orbitals energy at final time
    occupation_number::occupationType                           # Occupation number at final time
    Noccup::Vector{Int}                                         # Triplet (Nf,Np,Nv)
end

function makesolution(cache::RCACache, ::RCAMethod, solver::KohnShamSolver)
    occupation = make_occupation_number(solver.discretization, cache)
    RCASolution{typeof(cache.D), typeof(cache.U), typeof(cache.ϵ), typeof(occupation)}(cache.D, cache.U, cache.ϵ, occupation, cache.Noccup)
end


# Make occupation number

function make_occupation_number(::LDADiscretization, cache::RCACache)
    @unpack ϵ, n = cache
    index = findall(x->x ≠ 0, n)
    index_sort = sortperm(ϵ[index])
    new_index = index[index_sort]
    return [(string(i[2]+ i[1] -1, L_QUANTUM_LABELS[i[1]]), ϵ[i], n[i]) for i ∈ new_index]
end

function make_occupation_number(::LSDADiscretization, cache::RCACache)
    @unpack ϵ, n = cache
    index = findall(x->x ≠ 0, n)
    index_sort = sortperm(ϵ[index])
    new_index = index[index_sort]
    return [(string(i[2]+ i[1] -1, L_QUANTUM_LABELS[i[1]],SPIN_LABELS[i[3]]), ϵ[i], n[i]) for i ∈ new_index]
end

#####################################################################
#                          RCA STEPS
#####################################################################


function stopping_criteria!(cache::RCACache, ::RCAMethod, solver::KohnShamSolver)
    norm(cache.D - cache.Dprev) + abs(solver.energies[:Etot] - cache.energies_prev[:Etot])
end

function loopheader!(cache::RCACache, ::RCAMethod, solver::KohnShamSolver)
    @unpack energies = solver
    @unpack energies_prev, Dprev, D = cache
    cache.Dprev            .= cache.D
    cache.flag_degen        = false
    energies_prev[:Etot]    = energies[:Etot]
    energies_prev[:Ekin]    = energies[:Ekin]
    energies_prev[:Ecou]    = energies[:Ecou]
    energies_prev[:Ehar]    = energies[:Ehar]
end


function performstep!(cache::RCACache, method::RCAMethod, solver::KohnShamSolver)

    @unpack discretization, model, opts, energies = solver
    @unpack D, Dprev, U, ϵ, n = cache
    
    # STEP 1 : PREPARE THE EIGENVALUE PROBLEM
    prepare_eigenvalue_problem!(discretization, model, Dprev, opts.hartree)

    # STEP 2 : FIND ORBITALS AND CORRESPONFING ENERGIES
    find_orbital!(discretization, U, ϵ)

    # STEP 3 : FILL THE OCCUPATION NUMBER MATRIX ACCORDINGLY WITH THE AUFBAU PRINCIPLE
    #          The normalization of eigenvectors is made during this step to only 
    #          normalize the eigenvectors we need.
    aufbau!(cache, solver)

    if !cache.flag_degen
        # STEP 4 : COMPUTE A GUESS DENSITY
        density!(discretization, U, n, D)
        # STEP 5 : COMPUTE ALL ENERGIES
        energies[:Etot] = compute_total_energy(discretization, model, D, n, ϵ)
        energies[:Ekin] = compute_kinetic_energy(discretization, U, n)
        energies[:Ecou] = compute_coulomb_energy(discretization, U, n)
        energies[:Ehar] = compute_hartree_energy(discretization, D)
        !isthereExchangeCorrelation(model) || 
                (energies[:Eexc] = compute_exchangecorrelation_energy(discretization, model, D))
        #!(typeof(discretization) <: LSDADiscretization) || 
        #        (energies[:Ekincorr] = compute_kinetic_correlation_energy!(discretization, model, D))
    end

    # STEP  6 : COMPUTE THE NEW DENSITY
    update_density!(cache, method, solver)
end


function loopfooter!(::RCACache, ::RCAMethod, ::KohnShamSolver) end

function monitor(cache::RCACache, ::RCAMethod, ::KohnShamSolver)
    println("Relaxed Parameter : $(cache.t)")
    println("degeneracy ? : $(cache.flag_degen)")
    if cache.flag_degen
        println("Interpolation parameters : $(cache.tdegen)")
    end
end

#####################################################################
#                   CONSTANT DAMPLING ALGORITHM
#####################################################################

struct CDA{T} <: RCAMethod
    t::T
    function CDA(t::Real)
        @assert 0 ≤ t ≤ 1
        new{typeof(t)}(t)
    end
end

name(::CDA) = "CDA" 

function update_density!(cache::RCACache, ::CDA, solver::KohnShamSolver)
    @unpack t, D, Dprev, energies_prev = cache
    @unpack energies, discretization, model = solver

    if solver. niter > 0
        # UPDATE THE DENSITY
        @. D = t * D + (1 - t) * Dprev

        # UPDATE THE ENERGIES
        energies[:Etot] = t*energies[:Etot] + (1-t)*energies_prev[:Etot]
        energies[:Ekin] = t*energies[:Ekin] + (1-t)*energies_prev[:Ekin]
        energies[:Ecou] = t*energies[:Ecou] + (1-t)*energies_prev[:Ecou]
        energies[:Ehar] = compute_hartree_energy(discretization, D)
        if isthereExchangeCorrelation(model)
            energies[:Eexc] = compute_exchangecorrelation_energy(discretization, model, D)
        end
    end
    nothing
end




#####################################################################
#                   OPTIMAL DAMPLING ALGORITHM
#####################################################################

struct ODA{T<:Real} <: RCAMethod 
    t::T
    iter::Int
    function ODA(t::Real, iter::Int = 1)
        new{typeof(t)}(t,iter)
    end
end

name(::ODA) = "ODA" 

function update_density!(cache::RCACache, m::ODA, solver::KohnShamSolver)

    @unpack D, Dprev, energies_prev= cache
    @unpack discretization, model, energies = solver

    if solver.niter > 0
        
        if solver.niter < m.iter
            D .= cache.t * D + (1 - cache.t) * Dprev

            # UPDATE THE ENERGIES
            energies[:Etot] = cache.t*energies[:Etot] + (1-cache.t)*energies_prev[:Etot]
            energies[:Ekin] = cache.t*energies[:Ekin] + (1-tcache.t)*energies_prev[:Ekin]
            energies[:Ecou] = cache.t*energies[:Ecou] + (1-cache.t)*energies_prev[:Ecou]
            energies[:Ehar] = compute_hartree_energy(discretization, D)
            if isthereExchangeCorrelation(model)
                energies[:Eexc] = compute_exchangecorrelation_energy(discretization, model, D)
            end
            return nothing
        end
        
        energy_har01 = compute_hartree_mix_energy(discretization, D, Dprev)
        energy_har10 = compute_hartree_mix_energy(discretization, Dprev, D)

        # FIND THE OPTIMUM OCCUPATION
        cache.t, energies[:Etot] = find_minima_oda( energies[:Ekin], energies_prev[:Ekin], 
                                                    energies[:Ecou], energies_prev[:Ecou], 
                                                    energies[:Ehar], energies_prev[:Ehar], 
                                                    energy_har01, energy_har10,
                                                    D, Dprev, model, discretization)

        # UPDATE THE DENSITY
        D .= cache.t * D + (1 - cache.t) * Dprev

        # UPDATE THE ENERGIES
        t = cache.t
        energies[:Ekin] = t*energies[:Ekin] + (1-t)*energies_prev[:Ekin]
        energies[:Ecou] = t*energies[:Ecou] + (1-t)*energies_prev[:Ecou]
        energies[:Ehar] = t^2*energies[:Ehar] + (1-t)^2*energies_prev[:Ehar] + t*(1-t) * (energy_har01 + energy_har10)
        if isthereExchangeCorrelation(model)
            energies[:Eexc] = compute_exchangecorrelation_energy(discretization, model, D)
        end

    end
    nothing
end