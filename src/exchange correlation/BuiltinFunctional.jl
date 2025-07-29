"""
    abstract type BuiltinFunctional end

Exchange or Correlation functional implemented directly in AtomicKohnSham.jl. This has been
implemented to be compared with Libxc.jl and to use higher precision than Float64.
"""
abstract type BuiltinFunctional end

function BuiltinFunctional(identifier::Symbol; n_spin::Int = 1)
    n_spin in (1, 2) || error("n_spin needs to be 1 or 2")
    if identifier == :lda_x
        return SlaterXα(n_spin)
    else
        error("$identifier is not known for BuildinFunctional.")
    end
end

is_lda(f::BuiltinFunctional) = f.family == :lda

"""
    NoFunctional
"""
struct NoFunctional <: BuiltinFunctional
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
        zk::OptArray = C_NULL,
        vrho::OptArray = C_NULL)
    if func isa NoFunctional
        if !(zk == C_NULL)
            fill!(zk, 0)
        end

        if !(vrho == C_NULL)
            fill!(vrho, 0)
        end
        return
    end

    rho = ndims(rho) == 1 ? reshape(rho, :, 1) : rho

    if !(zk == C_NULL)
        fill!(zk, 0)
        if func.n_spin == 1
            for i in eachindex(rho)
                @views ρi = rho[i]
                zk[i] = eval_zk(func, ρi)
            end
        else
            for i in axes(rho, 2)
                @views ρi = rho[:, i]
                zk[i] = eval_zk(func, ρi...)
            end
        end
    end

    if !(vrho == C_NULL)
        fill!(vrho, 0)
        if func.n_spin == 1
            for i in eachindex(rho)
                @views ρi = rho[i]
                vrho[i] = eval_vrho(func, ρi)
            end
        else
            for i in axes(rho, 2)
                @views ρi = rho[:, i]
                vrho[1, i] = eval_vrho_up(func, ρi...)
                vrho[2, i] = eval_vrho_down(func, ρi...)
            end
        end
    end
end
