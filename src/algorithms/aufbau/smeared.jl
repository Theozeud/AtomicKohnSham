"""
    SmearedAufbau(; Temp)

Finite-temperature Aufbau occupation scheme using Fermi–Dirac smearing.

Occupation numbers are set to `n_i = g_i * f(ε_i; μ, Temp)`, where `g_i` is the
orbital's maximal degeneracy (`4l+2` or `2l+1`, see [`degeneracy`](@ref)),
`f(ε; μ, T) = 1/(1+exp((ε-μ)/T))` is the Fermi–Dirac distribution, and the
chemical potential `μ` is found by bisection so that `∑ᵢ n_i` matches the
model's electron count `N`.

`Temp` is the electronic smearing width in Hartree (energy units, not degrees).
Smaller `Temp` gives sharper, more Aufbau-like filling; `Temp` must be strictly
positive (the T=0 limit is exactly [`OptimizedAufbau`](@ref)).

!!! note
    This scheme only assigns fractional occupations from the Fermi–Dirac
    distribution; it does not add the electronic entropy term `-Temp*S` to the
    reported total energy. `Etot` is therefore the standard Kohn–Sham energy
    evaluated at those occupations, not the Mermin free energy.

# Parameters
- `Temp::Real`: Electronic smearing width (> 0), in Hartree.
"""
struct SmearedAufbau{T <: Real} <: Aufbau
    Temp::T
    function SmearedAufbau(; Temp::Real)
        Temp > 0 || error("Temp must be strictly positive.")
        new{typeof(Temp)}(Temp)
    end
end

"""
    SmearedAufbauCache

Cache structure for `SmearedAufbau`.

Stores the per-orbital maximal degeneracy and a reusable buffer for the
flattened orbital energies, to avoid repeated allocations during the SCF loop.
"""
struct SmearedAufbauCache{T <: Real} <: AufbauCache
    Temp::T
    degen::Vector{Int}
    ϵvec::Vector{T}
    postcomputations::Bool
    function SmearedAufbauCache(discretization::KSEDiscretization, model::KSEModel,
            aufbau::SmearedAufbau)
        ϵ = zero_orbitals_energies(discretization)
        degen = [degeneracy(discretization, i) for i in eachindex(ϵ)]
        ϵvec = zero(vec(ϵ))
        new{typeof(aufbau.Temp)}(aufbau.Temp, degen, ϵvec, false)
    end
end

function create_cache_aufbau(discretization::KSEDiscretization, model::KSEModel,
        aufbau::SmearedAufbau)
    SmearedAufbauCache(discretization, model, aufbau)
end

"""
    aufbau!(n, ϵ, U, model, discretization, aufbau::SmearedAufbauCache, niter; verbose=0)

Assign Kohn–Sham orbital occupations using Fermi–Dirac smearing.

The chemical potential `μ` is bisected so that the total Fermi–Dirac occupancy
matches `model.N`, then `n` is overwritten in place with `g_i*f(ε_i;μ,Temp)`.
"""
function aufbau!(n::AbstractArray{<:Real}, ϵ::AbstractArray{<:Real}, U::AbstractArray{<:Real},
        model::KSEModel, discretization::KSEDiscretization, aufbau::SmearedAufbauCache,
        niter::Int; verbose::UInt8 = 0x00)
    @unpack Temp, degen, ϵvec = aufbau
    copyto!(ϵvec, ϵ)
    N = model.N

    # Bracket the chemical potential between the lowest and highest orbital
    # energy, widened by a margin (in units of Temp) so N is attainable even
    # when it would sit right at the edge of the occupied spectrum.
    μlo = minimum(ϵvec) - 50*Temp
    μhi = maximum(ϵvec) + 50*Temp

    # Bisection: total occupancy is monotonic increasing in μ, so this always
    # converges regardless of the shape of the orbital-energy spectrum.
    for _ in 1:100
        μmid = (μlo + μhi)/2
        occ = zero(N)
        for i in eachindex(ϵvec)
            occ += degen[i]/(1 + exp((ϵvec[i] - μmid)/Temp))
        end
        if occ < N
            μlo = μmid
        else
            μhi = μmid
        end
    end
    μ = (μlo + μhi)/2

    for i in eachindex(ϵvec)
        n[i] = degen[i]/(1 + exp((ϵvec[i] - μ)/Temp))
    end
    nothing
end
