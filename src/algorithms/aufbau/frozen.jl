"""
    FrozenAufbau(n::Dict{String,<:Real})

Frozen occupation scheme (user-prescribed occupations).

`FrozenAufbau` enforces fixed Kohn–Sham orbital occupations specified by the user
through a dictionary mapping shell labels to occupation numbers. During SCF,
the occupations are copied verbatim at each iteration, ignoring the orbital
energies `ϵ` and the Aufbau ordering.

Keys must be shell strings , e.g. `"1s"`, `"2p"`, and, for spin-polarized discretizations,
include an explicit spin suffix such as `"1sUP"`/`"1sDoWN"`.

# Parameters
- `n::Dict{String,<:Real}`: Mapping `shell_string => occupation`.
  Examples: `Dict("1s" => 2.0, "2s" => 2.0, "2p" => 6.0)`.
"""
struct FrozenAufbau{T<:Real} <: Aufbau
    n::Dict{String,T}
    function FrozenAufbau(n::Dict{String,<:Real})
        T = promote_type(typeof.(values(n))...)
        nT = convert(Dict{String,T},n)
        new{T}(nT)
    end
end


"""
    FrozenAufbauCache

Cache structure for `FrozenAufbau`.

Stores the fixed occupation array `nfix` (in the discretization indexing
convention) built from the user-provided dictionary in `FrozenAufbau`.
"""
struct FrozenAufbauCache{N} <: AufbauCache
    nfix::N
    postcomputations::Bool
    function FrozenAufbauCache(discretization::KSEDiscretization, model::KSEModel,
                              aufbau::FrozenAufbau)
        postcomputations = false
        ndict = aufbau.n
        nfix = zero_occupation_numbers(discretization)
        layers_st = keys(ndict)
        layers_qn = parse_shell.(layers_st)
        if max(length.(layers_qn)...) < 3 && discretization.nspin > 1
            error("You have a spinned-model. Please provide occupation numbers with the
            spin index. For example \"1sUP\".")
        elseif max(length.(layers_qn)...) == 3 && discretization.nspin == 1
            error("You have a no spinned-model. Please provide occupation numbers without
            the spin index. For example \"1s\".")
        end
        for k ∈ layers_st
            lqn = parse_shell(k)
            n = lqn[1]
            l = lqn[2]
            a = if length(lqn) == 2
                (l+1, n-l)
            else
                (l+1, n-l, lqn[3])
            end
            nfix[a...] = ndict[k]
        end
        new{typeof(nfix)}(nfix, postcomputations)
    end
end

function create_cache_aufbau(discretization::KSEDiscretization, model::KSEModel,
        aufbau::FrozenAufbau)
    FrozenAufbauCache(discretization, model, aufbau)
end


"""
    aufbau!(n, ϵ, U, model, discretization, aufbau, niter; verbose=0)

Assign Kohn–Sham orbital occupations using a frozen (user-prescribed) scheme.

The occupation array `n` is overwritten in place by copying the cached fixed
occupations `aufbau.nfix`.
"""
function aufbau!(n::AbstractArray{<:Real}, ϵ::AbstractArray{<:Real}, U::AbstractArray{<:Real},
                model::KSEModel, discretization::KSEDiscretization, aufbau::FrozenAufbauCache,
                niter::Int; verbose::UInt8 = 0x00)
    copyto!(n, aufbau.nfix)
end
