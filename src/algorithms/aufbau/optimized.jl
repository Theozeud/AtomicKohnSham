"""
    OptimizedAufbau(; max_degen=2, tol=1e-3, handle_degen_spin_1iter=true)

Aufbau occupation scheme with explicit treatment of low-dimensional degeneracies.

This scheme fills orbitals following the Aufbau principle and
handles degenerate levels by solving a small variational problem
to minimize the total energy when the last shell is only partially filled.

# Parameters
- `max_degen::Int`: Maximum degeneracy dimension handled explicitly (default: 2).
- `tol::Real`: Energy tolerance to detect degeneracies.
- `handle_degen_spin_1iter::Bool`:
    If `true`, spin-only degeneracies (same spatial orbital, opposite spins)
    are treated symmetrically at the first SCF iteration.

"""
struct OptimizedAufbau{T<:Real} <: Aufbau
    max_degen::Int                  # Maximal Dimension of degeneracy considered
    tol::T                          # Tolerance for degenescence of orbitals energy
    handle_degen_spin_1iter::Bool   # True if at the first iteration, in case of degeneracy
                                    # between (nlσ₁) and (nlσ₂), this is handle differently.
    function OptimizedAufbau(;max_degen::Int = 2, tol::Real = 1e-3,
        handle_degen_spin_1iter::Bool = true)
        new{typeof(tol)}(max_degen, tol, handle_degen_spin_1iter)
    end
end


"""
    OptimizedAufbauCache

Cache structure for `OptimizedAufbau`.

Stores preallocated buffers used during the occupation procedure
to avoid repeated allocations and improve performance.

"""
mutable struct OptimizedAufbauCache{T<:Real, D <: AbstractArray{<:Real}, TF} <: AufbauCache
    max_degen::Int
    tol::T
    handle_degen_spin_1iter::Bool
    ϵvec::Vector{T}
    indices_sort::Vector{Int}
    indices_block::Vector{Int}
    degen_block::Vector{Int}
    D1 ::D
    D2::D
    Dbuf::D
    F::TF
    postcomputations::Bool
    energies::Energies{T}
    function OptimizedAufbauCache(discretization::KSEDiscretization, model::KSEModel,
        aufbau::OptimizedAufbau)
        @unpack max_degen, tol, handle_degen_spin_1iter = aufbau
        ϵ = zero_orbitals_energies(discretization)
        ϵvec = zero(vec(ϵ))
        indices_sort = zeros(Int,length(ϵ))
        indices_block = zeros(Int,max_degen)
        degen_block = zeros(Int,max_degen)

        dT = eltype(discretization)

        D1 = zero_density(discretization)
        D2 = zero_density(discretization)
        Dbuf = zero_density(discretization)
        F = if has_exchcorr(model)
            t -> begin
                @. Dbuf = t*D1 + (1-t)*D2
                compute_exc_energy(discretization, model, Dbuf)
            end
        else
            zero(dT)
        end

        postcomputations = false
        energies = Energies(dT)
        new{dT, typeof(D1), typeof(F)}(
            max_degen, tol, handle_degen_spin_1iter, ϵvec, indices_sort,
            indices_block, degen_block, D1, D2, Dbuf, F, postcomputations, energies)
    end
end

function create_cache_aufbau(discretization::KSEDiscretization, model::KSEModel,
        aufbau::OptimizedAufbau)
    OptimizedAufbauCache(discretization, model, aufbau)
end

"""
    aufbau!(n, ϵ, U, model, discretization, aufbau, niter; verbose=0)

Assign Kohn–Sham orbital occupations according to the Aufbau principle.

Orbitals are filled in increasing energy order, with automatic handling of
simple degeneracies (up to twofold). Fractional occupations are assigned when
the last occupied shell is partially filled. The occupation array `n` is
updated in place.
"""
function aufbau!(n::AbstractArray{<:Real}, ϵ::AbstractArray{<:Real},U::AbstractArray{<:Real},
    model::KSEModel, discretization::KSEDiscretization, aufbau::OptimizedAufbauCache,
    niter::Int; verbose::UInt8 = 0x00)
    @unpack ϵvec, indices_sort, indices_block, tol, max_degen, degen_block = aufbau
    copyto!(ϵvec,ϵ)
    sortperm!(indices_sort, ϵvec)
    fill!(n, 0)
    Nleft = model.N
    idx = 1
    nb_degen = 0
    aufbau.postcomputations = false
    while Nleft > 0 && idx ≤ length(indices_sort)
        # FIND ALL THE ORBITALS WITH THE SAME ENERGY
        nb_degen = collect_degenerate_block!(
            indices_block, indices_sort, ϵ, idx, tol, max_degen
        )
        idx += nb_degen
        # COMPUTE DEGENERACY
        total_degen = compute_total_degeneracy!(
            degen_block, discretization, indices_block, nb_degen
        )
        # ENOUGH ELECTRONS TO FILL ALL THE BLOCK
        if Nleft ≥ total_degen
            fill_full_block!(
                n, degen_block, indices_block, nb_degen
            )
            Nleft -= total_degen
        # NOT ENOUGH ELECTRONS FOR THE LAST SINGLE ORBITAL
        elseif nb_degen == 1
            fill_partial_single!(
                n, indices_block[1], Nleft
            )
            break
        # NOT ENOUGH ELECTRONS TO DISTRIBUTE BETWEEN TWO ORBITALS
        elseif nb_degen == 2
            resolve_twofold_degeneracy!(
                n, Nleft, U, aufbau, discretization, niter
            )
            break
        else
            error("You encouter a degeneracy $nb_degen. Degeneracy > 2 not implemented.")
        end
    end
    if verbose == 1
        println("Degeneracy of the last layer : $(nb_degen) between ")
        if nb_degen > 1
            for i ∈ 1:nb_degen
                idx = indices_block[i]
                if model.discretization.n_spin == 1
                    l,k,_ = convert_index(discretization, indices_block[1])
                    layer = shell_string((k+l,l-1))
                    println("layer : $layer")
                else
                    l,k,σ = convert_index(discretization, indices_block[1])
                    layer = shell_string((k+l,l-1),σ)
                    println("layer : $layer")
                end
            end
        end
    end
end


"""
Collect orbitals with identical energy (up to two).
Uses a cache buffer to avoid allocations.
"""
function collect_degenerate_block!(
    indices_block::Vector{Int}, indices_sort::Vector{Int}, ϵ::AbstractArray{<:Real},
    idx::Int, tol::Real, max_degen::Real
)
    indices_block[1] = indices_sort[idx]
    nb_degen = 1
    ϵixd = ϵ[indices_block[1]]
    j = idx+1
    i = 2
    while j ≤ length(indices_sort) &&
          abs(ϵ[indices_sort[j]] - ϵixd) < tol && nb_degen < max_degen
        indices_block[i] = indices_sort[j]
        j += 1
        i += 1
        nb_degen += 1
    end
    return nb_degen
end

"""
Compute total degeneracy of a block of orbitals.
"""
function compute_total_degeneracy!(
    degen_block::Vector{Int}, discretization::KSEDiscretization,
    indices_block::Vector{<:Real}, nb_degen::Int
)
    fill!(degen_block, 0)
    total_degen = 0
    for i ∈ 1:nb_degen
        idx = indices_block[i]
        degen_block[i] = degeneracy(discretization, idx)
        total_degen += degen_block[i]
    end
    return total_degen
end

"""
Fill completely all orbitals of a degenerate block.
"""
function fill_full_block!(
    n::AbstractArray{<:Real}, degen_block::Vector{Int}, indices_block::Vector{Int},
    nb_degen::Int
)
    for i ∈ 1:nb_degen
        idx = indices_block[i]
        n[idx] = degen_block[i]
    end
end

"""
Partially fill a single orbital.
"""
function fill_partial_single!(
    n::AbstractArray{<:Real}, idx::Int, Nleft::Real
)
    n[idx] = Nleft
end

"""
Resolve a twofold degeneracy by energy minimization.
Includes a special symmetric spin case.
"""
function resolve_twofold_degeneracy!(
    n::AbstractArray{<:Real}, Nleft::Int, U::AbstractArray{<:Real},
    aufbau::OptimizedAufbauCache, discretization::KSEDiscretization,
    niter::Int
)
    @unpack indices_block, handle_degen_spin_1iter, degen_block, energies,
            D1, D2, F = aufbau

    l1,k1,σ1 = convert_index(discretization, indices_block[1])
    l2,k2,σ2 = convert_index(discretization, indices_block[2])
    if handle_degen_spin_1iter && niter == 0 && l1==l2 && k1==k2 && σ1≠σ2
        n[indices_block[1]] = Nleft/2
        n[indices_block[2]] = Nleft/2
        return
    end

    aufbau.postcomputations = true

    # COMPUTE ENERGIES FOR ONE EXTREMA
    degen_first = degen_block[1]
    if degen_first ≥ Nleft
        n[indices_block[1]] = Nleft
        n[indices_block[2]] = zero(Nleft)
        n1_0 = Nleft
        n2_0 = zero(Nleft)
    else
        n[indices_block[1]] = degen_first
        n[indices_block[2]] = Nleft - degen_first
        n1_0 = degen_first
        n2_0 = Nleft - degen_first
    end
    density!(discretization, U, n, D1)
    energy_kin0 = compute_kinetic_energy(discretization, U, n)
    energy_cou0 = compute_coulomb_energy(discretization, U, n)
    energy_har0 = compute_hartree_energy(discretization, D1)

    # COMPUTE ENERGIES FOR THE OTHER EXTREMA
    degen_second = degen_block[2]
    if degen_second ≥ Nleft
        n[indices_block[1]] = zero(Nleft)
        n[indices_block[2]] = Nleft
        n1_1 = zero(Nleft)
        n2_1 = Nleft
    else
        n[indices_block[1]] = Nleft - degen_second
        n[indices_block[2]] = degen_second
        n1_1 = Nleft - degen_second
        n2_1 = degen_second
    end
    density!(discretization, U, n, D2)
    energy_kin1 = compute_kinetic_energy(discretization, U, n)
    energy_cou1 = compute_coulomb_energy(discretization, U, n)
    energy_har1 = compute_hartree_energy(discretization, D2)

    # COMPUTE THE QUADRATIC TERM RELATED TO HARTREE ENERGY
    energy_har01 = compute_hartree_mix_energy(discretization, D1, D2)

    # FIND THE OPTIMUM OCCUPATION
    tdegen, energies.Etot = line_search_energy(
        energy_kin0, energy_kin1,
        energy_cou0, energy_cou1,
        energy_har0, energy_har1,
        energy_har01, F)

    # UPDATE THE OCCUPATION NUMBERS
    n[indices_block[1]] = tdegen * n1_0 + (1-tdegen) * n1_1
    n[indices_block[2]] = tdegen * n2_0 + (1-tdegen) * n2_1

    # UPDATE THE DENSITY
    @. D1 = tdegen * D1 + (1 - tdegen) * D2

    # UPDATE THE ENERGIES
    t = tdegen
    energies.Ekin = t*energy_kin0 + (1-t)*energy_kin1
    energies.Ecou = t*energy_cou0 + (1-t)*energy_cou1
    energies.Ehar = t^2*energy_har0 + (1-t)^2*energy_har1 +
                        2 * t * (1-t) * energy_har01
    energies.Eexc = energies.Etot-energies.Ekin-energies.Ecou-energies.Ehar
end
