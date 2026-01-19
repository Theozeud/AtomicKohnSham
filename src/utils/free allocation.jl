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
