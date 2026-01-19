antiadjoint(z::Complex) = -real(z) + imag(z)
antiadjoint(z::Real) = -z

function _mul!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
        D::AbstractMatrix, tmp::AbstractMatrix)
    # COMPUTE B*C*D -> A USING tmp AS TEMPORY VARIABLES
    mul!(tmp, C, D)
    mul!(A, B, tmp)
end

function _commutator!(
        A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, tmp::AbstractMatrix)
    # COMPUTE [B,C] -> A USING tmp AS TEMPORY VARIABLES
    mul!(A, B, C)
    mul!(tmp, C, B)
    @. A -= tmp
end


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
