# ===================================================================
#                             Density
# ===================================================================
"""
    density!(discretization, U, n, D)

Assemble the one-particle density matrix `D` from Kohn–Sham orbitals `U`
and occupation numbers `n`.

The density matrix is constructed as
    D = ∑ₗₖσ nₗₖσ |Uₗₖσ⟩⟨Uₗₖσ|,
where the sum runs over angular momentum, radial index, and spin.
The result is symmetric and stored in-place in `D`.
"""
function density!(discretization::KSEDiscretization{T},
                 U::AbstractArray{<:Real},
                 n::AbstractArray{<:Real},
                 D::AbstractArray{<:Real}) where T
    @unpack lₕ, nₕ, Nₕ, nspin = discretization
    fill!(D, zero(T))
    @inbounds for σ in 1:nspin
        @inbounds for k in 1:nₕ
            @inbounds for l in 1:(lₕ + 1)
                if !iszero(n[l, k, σ])
                    @inbounds for i in 1:Nₕ
                        val = n[l, k, σ] * U[i, k, l, σ]
                        @inbounds @simd for j in 1:i
                            D[i, j, σ] += val * U[j, k, l, σ]
                        end
                    end
                end
            end
        end
    end
    @inbounds for i in 1:Nₕ
        @inbounds @simd for j in 1:(i - 1)
            D[j, i, :] .= D[i, j, :]
        end
    end
    nothing
end

# ===================================================================
#                       Evaluation of Density
# ===================================================================
"""
    density!(discretization, U, n, D)

Assemble the one-particle density matrix `D` from Kohn–Sham orbitals `U`
and occupation numbers `n`.

The density matrix is constructed as
    D = ∑ₗₖσ nₗₖσ |Uₗₖσ⟩⟨Uₗₖσ|,
where the sum runs over angular momentum, radial index, and spin.
The result is symmetric and stored in-place in `D`.
"""
function eval_density(discretization::KSEDiscretization,
                     D::AbstractMatrix{<:Real},
                     x::Real)
    @unpack basis, cache = discretization
    @unpack buf1, buf2 = cache.evalw
    localisation_x = findindex(basis.mesh, x)
    I = basis.cells_to_indices[localisation_x]
    @views eval_basis = buf2[I]
    evaluate!(eval_basis, basis, I, x)
    @views tv = buf1[I]
    @views Dview = D[I, I]
    mul!(tv, Dview, eval_basis)
    return 1/(4π*x^2) * dot(eval_basis, tv)
end

"""
    eval_density(discretization, D, X)

Evaluate the radial electronic density at all positions in `X`
from the density matrix `D`.

Returns a vector containing ρ(x) for each x ∈ X.
"""
function eval_density(discretization::KSEDiscretization,
                     D::AbstractMatrix{<:Real},
                     X::AbstractVector{<:Real})
    newT = promote_type(eltype(D), eltype(X))
    ρ = zeros(newT, length(X))
    eval_density!(ρ, discretization, D, X)
    ρ
end

"""
    eval_density!(ρ, discretization, D, X)

Evaluate the radial electronic density at positions `X` from the density
matrix `D` and store the result in `ρ`.

This in-place version avoids allocations and reuses internal work buffers.
"""
function eval_density!(ρ::AbstractVector{<:Real},
                      discretization::KSEDiscretization,
                      D::AbstractMatrix{<:Real},
                      X::AbstractVector{<:Real})
    @unpack basis, cache = discretization
    @unpack buf1, buf2 = cache.evalw
    cache_Pϕx = _cache_Pϕx(basis, first(X))
    @inbounds for k in eachindex(X)
        xk = X[k]
        localisation_xk = findindex(basis.mesh, xk)
        Ik = basis.cells_to_indices[localisation_xk]
        @views eval_basis = buf2[Ik]
        evaluate!(eval_basis, basis, Ik, xk, cache_Pϕx)
        @views tv = buf1[Ik]
        @views Dk = D[Ik, Ik]
        mul!(tv, Dk, eval_basis)
        ρ[k] = 1/(4π*xk^2) * dot(eval_basis, tv)
    end
    ρ
end
