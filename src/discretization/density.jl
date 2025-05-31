#--------------------------------------------------------------------
#                             Density
#--------------------------------------------------------------------
function density!(  discretization::LDADiscretization, 
                    U::AbstractArray{<:Real}, 
                    n::AbstractMatrix{<:Real}, 
                    D::AbstractMatrix{<:Real})
    @unpack lₕ, nₕ, Nₕ, elT  = discretization
    fill!(D, zero(elT))
    @inbounds for k ∈ 1:nₕ
        @inbounds for l ∈ 1:lₕ+1   
            if !iszero(n[l,k])
                @inbounds for i ∈ 1:Nₕ
                    val = n[l,k] * U[i,k,l] 
                    @inbounds @simd for j ∈ 1:i
                        D[i,j] += val * U[j,k,l]
                    end
                end
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


function density!(  discretization::LDADiscretization, 
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


function compute_density(discretization::LDADiscretization, D::AbstractMatrix{<:Real}, x::Real)
    @unpack basis, cache = discretization
    @unpack tmp_vect, tmp_C = cache
    localisation_x = findindex(basis.mesh, x)
    I = basis.cells_to_indices[localisation_x]
    @views eval_basis = tmp_C[I]
    evaluate!(eval_basis, basis.)
    #=
    @inbounds for (n,i) ∈ enumerate(I)
    eval_basis[n] = basis(i,x)
    end
    =#
    @views tv = tmp_vect[I]
    @views Dview = D[I,I]
    mul!(tv,Dview,eval_basis)
    return 1/(4π*x^2) * dot(eval_basis,tv)
end


#--------------------------------------------------------------------
#                          Density Matrix
#--------------------------------------------------------------------

function density_matrix!(   discretization::LDADiscretization, 
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
