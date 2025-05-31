function _spzeros(T::Type, n::Int, m::Int, p::Int)
    A = Vector{SparseMatrixCSC{T, Int}}(undef, p)
    for i in 1:p
        A[i] = spzeros(T, n, m)
    end
    A
end

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