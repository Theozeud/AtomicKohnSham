"""
    abstract type BuiltinFunctional end

Exchange or Correlation functional implemented directly in AtomicKohnSham.jl. This has been
implemented to be compared with Libxc.jl and to use higher precision than Float64.
"""
abstract type BuiltinFunctional end

# Unlike Libxc.OptArray (Union{Ptr{Nothing}, Array{Float64}} -- Libxc calls
# into a Float64-only C library), BuiltinFunctional is pure Julia and must
# accept any real element type to honor the "higher precision than Float64"
# claim above.
const OptRealArray = Union{Ptr{Nothing}, AbstractArray{<:Real}}

function BuiltinFunctional(identifier::Symbol;nspin::Int = 1)
   nspin in (1, 2) || error("n_spin needs to be 1 or 2")
    if identifier == :lda_x
        return SlaterXα(nspin)
    elseif identifier == :lda_c_pw
        return PerdewWang(nspin)
    else
        error("$identifier is not known for BuildinFunctional.")
    end
end

is_lda(f::BuiltinFunctional) = f.family == :lda

# Every BuiltinFunctional implemented in AtomicKohnSham.jl today (SlaterXα,
# PerdewWang, NoFunctional) is LDA-only -- GGA builtins are a separate,
# deferred task (see dev/gga/PLAN.md). NoFunctional has no `family` field, so
# this can't reuse `is_lda`'s `f.family == :lda` pattern.
is_gga(f::BuiltinFunctional) = false

"""
    NoFunctional
"""
@Base.kwdef struct NoFunctional <: BuiltinFunctional
   n_spin::Int
end

function evaluate_functional(
        func::BuiltinFunctional;
        rho::AbstractArray{<:Real},
        derivatives = 0:1)
    if ndims(rho) > 1 && size(rho, 1) != func.n_spin
        error("First axis for multidimensional rho array should be equal " *
              "to the number of spin components (== $(func.n_spin))")
    end

    if derivatives == 0
        zk = zeros(eltype(rho), length(rho))
        evaluate_functional!(func; rho = rho, zk = zk, vrho = C_NULL)
        return (zk = zk,)
    elseif derivatives == 1
        vrho = if func.n_spin == 1
            zeros(eltype(rho), length(rho))
        else
            zeros(eltype(rho), 2, length(rho))
        end
        evaluate_functional!(func; rho = rho, zk = C_NULL, vrho = vrho)
        return (vhro = vrho,)
    elseif derivatives == 0:1
        zk = zeros(eltype(rho), length(rho))
        vrho = if func.n_spin == 1
            zeros(eltype(rho), length(rho))
        else
            zeros(eltype(rho), 2, length(rho))
        end
        evaluate_functional!(func; rho = rho, zk = zk, vrho = vrho)
        return (vhro = vrho, zk = zk)
    end
end

function evaluate_functional!(
        func::BuiltinFunctional;
        rho::AbstractArray{<:Real},
        zk::OptRealArray = C_NULL,
        vrho::OptRealArray = C_NULL,
        sigma = nothing,
        vsigma::OptRealArray = C_NULL)
    # sigma/vsigma: accepted (and ignored beyond zeroing vsigma) so that a
    # BuiltinFunctional -- always LDA-only here -- can stand in for the
    # non-GGA half of a GGA model (e.g. NoFunctional correlation alongside a
    # GGA exchange); see has_gga/KSEModel's rung check in physics/models.jl.
    if func isa NoFunctional
        if !(zk == C_NULL)
            fill!(zk, 0)
        end

        if !(vrho == C_NULL)
            fill!(vrho, 0)
        end

        if !(vsigma == C_NULL)
            fill!(vsigma, 0)
        end
        return
    end

    rho = ndims(rho) == 1 ? reshape(rho, :, 1) : rho

    # Density should be non-negative by construction, but a FEM/quadrature
    # approximation of ρ can dip slightly (or, for a badly-conditioned
    # discretization -- too few, too-high-order elements over too large a
    # domain -- wildly) negative. eval_zk/eval_vrho below raise ρ to
    # non-integer powers (e.g. ρ^(1/3)), which throws a DomainError for any
    # negative real input rather than degrading gracefully. Clamping at the
    # physical floor of zero (a functional's only sane value in vacuum
    # anyway, see the iszero(ρ) guards elsewhere in this file/PerdewWang.jl)
    # turns that hard crash into the same "no contribution here" behavior a
    # tiny negative artifact should have had in the first place.
    _nonneg(x) = max(x, zero(x))

    if !(zk == C_NULL)
        fill!(zk, 0)
        if func.n_spin == 1
            for i in eachindex(rho)
                zk[i] = eval_zk(func, _nonneg(rho[i]))
            end
        else
            for i in axes(rho, 2)
                @views ρi = rho[:, i]
                zk[i] = eval_zk(func, _nonneg(ρi[1]), _nonneg(ρi[2]))
            end
        end
    end

    if !(vrho == C_NULL)
        fill!(vrho, 0)
        if func.n_spin == 1
            for i in eachindex(rho)
                vrho[i] = eval_vrho(func, _nonneg(rho[i]))
            end
        else
            for i in axes(rho, 2)
                @views ρi = rho[:, i]
                vrho[1, i] = eval_vrho_up(func, _nonneg(ρi[1]), _nonneg(ρi[2]))
                vrho[2, i] = eval_vrho_down(func, _nonneg(ρi[1]), _nonneg(ρi[2]))
            end
        end
    end
end
