#########################################################################################
#                                     NUMERICAL TOOLS
#########################################################################################

# Solving linear problem
function solve_linear_problem(A::AbstractMatrix{<:Real}, b::AbstractVector{<:Real})
    return A\b
end

# Solve Generalized Eigenvalue problem
function solve_generalized_eigenvalue_problem(A::AbstractMatrix{<:Real}, B::AbstractMatrix{<:Real}, n::Int)
    #C = sqrt(inv(B))
    #λ, U = real.(eigs(A,B; which = :SR, nev = n))
    #λ, U
    C = sqrt(inv(B))
    λ, U = eigen(C*A*C)
    #λ, U = eigen(A,B)
    λ[1:n], C*U[:,1:n]
 
end

# Approximate integral
function approximate_integral(f, domain; method = QuadGKJL(), reltol = 1e-3, abstol = 1e-3)
    prob = IntegralProblem((x,p) -> f(x), domain)
    solve(prob, method; reltol = 1e-8, abstol = 1e-8).u
end


#########################################################################################
#                             FREE ALLOCATIONS FUNCTION
#########################################################################################

antiadjoint(z::Complex) = -real(z) + imag(z)
antiadjoint(z::Real) = -z

function _mul!( A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, 
                D::AbstractMatrix, tmp::AbstractMatrix)
    # COMPUTE B*C*D -> A USING tmp AS TEMPORY VARIABLES
    mul!(tmp,C,D)
    mul!(A,B,tmp)
end

function _commutator!( A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, tmp::AbstractMatrix)
    # COMPUTE [B,C] -> A USING tmp AS TEMPORY VARIABLES
    mul!(A,B,C)
    mul!(tmp,C,B)
    @. A -= tmp
end


function mul_block_by_block!(Γout::BlockDiagonal{<:Number}, M::AbstractMatrix, Γin::BlockDiagonal{<:Number})
    if nbblock(Γout) ≠ nbblock(Γin)
        throw(DimensionMismatch(lazy"Γout has $(nbblock(Γout)) blocks but Γin has $(nbblock(Γin))"))
    end
    for i ∈ 1:nbblock(Γout)
        @views Γin_i = blocks(Γin)[i]
        @views Γout_i = blocks(Γout)[i]
        mul!(Γout_i, M, Γin_i)
    end
    Γout
end


function mul_block_by_block!(Γout::BlockDiagonal{<:Number}, Γin::BlockDiagonal{<:Number}, M::AbstractMatrix,)
    if nbblock(Γout) ≠ nbblock(Γin)
        throw(DimensionMismatch(lazy"Γout has $(nbblock(Γout)) blocks but Γin has $(nbblock(Γin))"))
    end
    for i ∈ 1:nbblock(Γout)
        @views Γin_i = blocks(Γin)[i]
        @views Γout_i = blocks(Γout)[i]
        mul!(Γout_i, Γin_i, M)
    end
    Γout
end

function _copy_vec_to_mat!( M::AbstractMatrix, ir_dest::AbstractRange{Int}, jr_dest::AbstractRange{Int},
                            V::AbstractVector, ir_src::AbstractRange{Int})        
    if length(ir_dest) * length(jr_dest) != length(ir_src)
        throw(ArgumentError(LazyString("source and destination must have same size (got ",
                                length(ir_dest)*length(jr_dest)," and ",length(ir_dest),")")))
    end
    @boundscheck checkbounds(M, ir_dest, jr_dest)
    @boundscheck checkbounds(V, ir_src)
    isrc = first(ir_src)
    for jdest in jr_dest
        for idest in ir_dest
            M[idest,jdest] = V[isrc]
            isrc += step(ir_src)
        end
    end
    return M
end

function _copy_mat_to_vec!( V::AbstractVector, ir_dest::AbstractRange{Int},
                            M::AbstractMatrix, ir_src::AbstractRange{Int}, jr_src::AbstractRange{Int})        
    if length(ir_src) * length(jr_src) != length(ir_dest)
        throw(ArgumentError(LazyString("source and destination must have same size (got ",
                                length(ir_dest)," and ",length(ir_src)*length(jr_src),")")))
    end
    @boundscheck checkbounds(V, ir_dest)
    @boundscheck checkbounds(M, ir_src, jr_src)
    idest = first(ir_dest)
    for jsrc in jr_src
        for isrc in ir_src
            V[idest] = M[isrc,jsrc]
            idest += step(ir_dest)
        end
    end
    return M
end

function remove_trace!(A::AbstractMatrix, B::AbstractMatrix)
    A .= B
    trace = tr(B)/size(B,1)
    @threads for i ∈ axes(M,1)
        @inbounds A[i,i] -= trace
    end
    nothing
end


#########################################################################################
#                             SPARSE TENSOR TOOLS
#########################################################################################

function _spzeros(T::Type, n::Int, m::Int, p::Int)
    A = Vector{SparseMatrixCSC{T, Int}}(undef, p)
    for i in 1:p
        A[i] = spzeros(T, n, m)
    end
    A
end