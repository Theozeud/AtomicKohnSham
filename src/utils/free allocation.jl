function flexible_zeros(T::Type, dims::NTuple{N, Int}, lastdim::Int) where {N}
    if lastdim == 1
        return zeros(T, dims...)
    else
        return zeros(T, dims..., lastdim)
    end
end

function flexible_zeros(T::Type, firstdim::Int, dims::NTuple{N, Int}) where {N}
    if firstdim == 1
        return zeros(T, dims...)
    else
        return zeros(T, firstdim, dims...)
    end
end

function tensor_matrix_dict!(B::AbstractVector{<:Real},
        D::AbstractMatrix{<:Real},
        F::Dict{Tuple{Int, Int, Int}, <:Real})
    fill!(B, zero(eltype(B)))
    @inbounds for ((i, j, m), F_ijm) in F
        B[m] += D[i, j] * F_ijm
    end
    nothing
end

function tensor_matrix_dict!(B::AbstractVector{<:Real},
        DUP::AbstractMatrix{<:Real},
        DDOWN::AbstractMatrix{<:Real},
        F::Dict{Tuple{Int, Int, Int}, <:Real})
    fill!(B, zero(eltype(B)))
    @inbounds for ((i, j, m), F_ijm) in F
        B[m] += (DUP[i, j] + DDOWN[i, j]) * F_ijm
    end
    nothing
end

function tensor_vector_dict!(B::AbstractMatrix{<:Real}, D::AbstractVector{<:Real},
        F::Dict{Tuple{Int, Int, Int}, <:Real})
    fill!(B, zero(eltype(B)))
    @inbounds for ((i, j, m), F_ijm) in F
        B[i, j] += D[m] * F_ijm
    end
    nothing
end
