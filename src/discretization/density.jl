#--------------------------------------------------------------------
#                             Density
#--------------------------------------------------------------------
function density!(  discretization::KSEDiscretization,
                    U::AbstractArray{<:Real},
                    n::AbstractArray{<:Real},
                    D::AbstractArray{<:Real})
    @unpack lₕ, nₕ, Nₕ, n_spin  = discretization
    elT = eltype(discretization)
    fill!(D, zero(elT))
    @inbounds for k ∈ 1:nₕ
        @inbounds for σ ∈ 1:n_spin
            @inbounds for l ∈ 1:lₕ+1
                if !iszero(n[l,k,σ])
                    @inbounds for i ∈ 1:Nₕ
                        val = n[l,k] * U[i,k,l,σ]
                        @inbounds @simd for j ∈ 1:i
                            D[i,j,σ] += val * U[j,k,l,σ]
                        end
                    end
                end
            end
        end
    end
    @inbounds for i in 1:Nₕ
        @inbounds @simd for j in 1:i-1
            D[j,i,:] .= D[i,j,:]
        end
    end
    nothing
end


function density!(  discretization::KSEDiscretization,
                    Γ::BlockDiagonal{<:Real, <:AbstractMatrix{<:Real}},
                    D::AbstractMatrix{<:Real})
    @unpack lₕ, Nₕ, elT  = discretization
    fill!(D, zero(elT))
    @inbounds @simd for l ∈ 1:lₕ+1
        @views Γl = blocks(Γ)[l]
        @inbounds for i ∈ 1:Nₕ
            @inbounds for j ∈ 1:i
                D[i,j] += (2*l + 1)*Γl[i,j]
            end
        end
    end
    @inbounds for i in 1:Nₕ
        @inbounds @simd for j in 1:i-1
            D[j,i] = D[i,j]
        end
    end
    nothing
end

#--------------------------------------------------------------------
#                       Evaluation of Density
#--------------------------------------------------------------------
function eval_density(  discretization::KSEDiscretization,
                        D::AbstractMatrix{<:Real},
                        x::Real)
    @unpack basis, cache = discretization
    @unpack tmp_vect, tmp_C = cache
    localisation_x = findindex(basis.mesh, x)
    I = basis.cells_to_indices[localisation_x]
    @views eval_basis = tmp_C[I]
    evaluate!(eval_basis, basis, I, x)
    @views tv = tmp_vect[I]
    @views Dview = D[I,I]
    mul!(tv,Dview,eval_basis)
    return 1/(4π*x^2) * dot(eval_basis,tv)
end


function eval_density(  discretization::KSEDiscretization,
                        D::AbstractMatrix{<:Real},
                        X::AbstractVector{<:Real})
    newT = promote_type(eltype(D), eltype(X))
    ρ = zeros(newT, length(X))
    eval_density!(ρ, discretization, D, X)
    ρ
end


function eval_density!( ρ::AbstractVector{<:Real},
                        discretization::KSEDiscretization,
                        D::AbstractMatrix{<:Real},
                        X::AbstractVector{<:Real})
    @unpack basis, cache = discretization
    @unpack tmp_vect, tmp_C = cache
    cache_Pϕx = _cache_Pϕx(basis, first(X))
    @inbounds for k ∈ eachindex(X)
        xk = X[k]
        localisation_xk = findindex(basis.mesh, xk)
        Ik = basis.cells_to_indices[localisation_xk]
        @views eval_basis = tmp_C[Ik]
        evaluate!(eval_basis, basis, Ik, xk, cache_Pϕx)
        @views tv = tmp_vect[Ik]
        @views Dk = D[Ik,Ik]
        mul!(tv, Dk, eval_basis)
        ρ[k] = 1/(4π*xk^2) * dot(eval_basis,tv)
    end
    ρ
end


function optimized_eval_density!(ρ::AbstractVector{<:Real},
                                 discretization::KSEDiscretization,
                                 D::AbstractMatrix{<:Real},
                                 X::AbstractVector{<:Real})
    @unpack basis, fem_integration_method, mesh = discretization
    Qgenx = fem_integration_method.Qgenx
    fill!(ρ, 0)
    idxmesh = findindex(mesh, first(X))
    Ib = basis.cells_to_indices[idxmesh]
    Ig = basis.cells_to_generators[idxmesh]
    @views DIb = D[Ib, Ib]
    @views QgenxIg = Qgenx[Ig, Ig, :]
    @tensor ρ[k] = DIb[i,j] * QgenxIg[i,j,k]
    @.ρ /= X^2
    @. ρ *= 1/4π
    ρ
end

#--------------------------------------------------------------------
#                          Density Matrix
#--------------------------------------------------------------------

function density_matrix!(   discretization::KSEDiscretization,
                            U::AbstractArray{<:Real},
                            n::AbstractMatrix{<:Real},
                            Γ::BlockDiagonal{<:Real, <:AbstractMatrix{<:Real}})
    @unpack lₕ, nₕ, Nₕ, elT  = discretization
    @inbounds for l ∈ 1:lₕ+1
        @views Γl = blocks(Γ)[l]
        fill!(Γl, zero(elT))
        @inbounds for k ∈ 1:nₕ
            if !iszero(n[l,k])
                @inbounds for i ∈ 1:Nₕ
                    val = n[l,k]/(2*l-1) * U[l,i,k]
                    @inbounds @simd for j ∈ 1:i
                        Γl[i,j] += val * U[l,j,k]
                    end
                end
            end
        end
        @inbounds for i in 1:Nₕ
            @inbounds @simd for j in 1:i-1
                Γl[j,i] = Γl[i,j]
            end
        end
    end
    nothing
end
