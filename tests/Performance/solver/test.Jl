using KohnShamResolution: LDADiscretization, tensor_matrix_dict!, tensor_vector_dict!

using SparseArrays
function symmetrize_sparse!(A::SparseMatrixCSC{Float64,Int})
    @assert size(A, 1) == size(A, 2) "Matrix must be square"
    n = size(A, 1)
    for j in 1:n
        for idx in A.colptr[j]:(A.colptr[j+1]-1)
            i = A.rowval[idx]
            if i != j  # skip diagonal
                # Find symmetric entry (i, j) and (j, i)
                # locate position of (j, i)
                @views Av = A.rowval[A.colptr[i]:(A.colptr[i+1]-1)]
                k = findfirst(x -> x == j, Av)
                if k !== nothing
                    k = k + A.colptr[i] - 1
                    avg = (A.nzval[idx] + A.nzval[k]) / 2
                    A.nzval[idx] = avg
                    A.nzval[k] = avg
                end
            end
        end
    end
    return A
end

function scale_sparse!(A::SparseMatrixCSC{T}, coeff::T) where T
    @inbounds for i in eachindex(A.nzval)
        A.nzval[i] *= coeff
    end
    return A
end

function hartree_matrix2!(discretization::LDADiscretization, D::AbstractMatrix{<:Real}, coeff::Real = 1.0)
    @unpack Rmax, matrices, cache = discretization
    @unpack A, M₀, F, Hartree = matrices
    @unpack tmp_MV, tmp_B, tmp_C = cache
    tensor_matrix_dict!(tmp_B, D, F)

    @show @allocated tmp_C .= A\tmp_B
    @tensor newCrho = D[i,j] * M₀[i,j]
    tensor_vector_dict!(tmp_MV, tmp_C, F)
    @show @allocated @. Hartree = tmp_MV + newCrho/Rmax * M₀
    @show @allocated scale_sparse!(Hartree, coeff)
    @show @allocated symmetrize_sparse!(Hartree)

    nothing
end


hartree_matrix2!(solver.discretization, solver.cache.D)