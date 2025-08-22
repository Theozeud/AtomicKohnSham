using UnPack
using LinearAlgebra
using SparseArrays

struct PolySet{T, typeData <: AbstractVecOrMat}
    coeffs::typeData
    function PolySet(polys::AbstractVecOrMat)
        new{eltype(polys), typeof(polys)}(polys)
    end
end

function Base.show(io::IO, ::MIME"text/plain", ps::PolySet)
    n, degmax_plus1 = size(ps.coeffs)
    println(io, "PolySet with $n polynomials (degmax = $(degmax_plus1 - 1))")
    show(io, MIME"text/plain"(), ps.coeffs)
end

Base.eltype(p::PolySet) = eltype(p.coeffs)
Base.getindex(p::PolySet, i::Int, j::Int) = p.coeffs[i, j]

Base.size(ps::PolySet) = size(ps.coeffs)
Base.size(ps::PolySet, n::Int) = size(ps.coeffs, n)
npolys(ps::PolySet) = size(ps.coeffs, 1)
maxdeg(ps::PolySet) = size(ps.coeffs, 2) - 1

function allocate_polyset(T::DataType, nbpoly::Int, degmax::Int)
    PolySet(zeros(T, nbpoly, degmax+1))
end

function Base.getindex(ps::PolySet, row::Int)
    @views v = ps.coeffs[row, :]
    return PolySet(v')
end

function Base.getindex(ps::PolySet, rows::AbstractVector{Int})
    @views v = ps.coeffs[rows, :]
    return PolySet(v)
end

function SparseArrays.sparse(ps::PolySet)
    sparseps = sparse(ps.coeffs)
    PolySet(sparseps)
end

function evaluate(ps::PolySet, x::Number)
    NewT = promote_type(typeof(x), eltype(ps))
    y = zeros(NewT, size(ps, 1))
    evaluate!(y, ps, x)
end

function evaluate(ps::PolySet, x::AbstractVector)
    NewT = promote_type(eltype(x), eltype(ps))
    y = zeros(NewT, size(ps, 1), length(x))
    evaluate!(y, ps, x)
end

"""
    evaluate!(y::AbstractVector, ps::PolySet, x::Number)

Evaluates a set of polynomials represented by `ps` at the scalar point `x`, writing the result in-place into the vector `y`.

# Arguments
- `y`: Output vector (modified in-place), with length equal to the number of polynomials. After evaluation, `y[i]` contains the value of the i-th polynomial at `x`.
- `ps`: A `PolySet` object storing polynomial coefficients. Each row corresponds to one polynomial; columns represent increasing powers of `x`.
- `x`: The scalar value at which all polynomials are evaluated.

# Notes
- Uses Horner's method for efficient polynomial evaluation.
- No heap allocation is performed during evaluation.
"""
function evaluate!(y::AbstractVector, ps::PolySet, x::Number)
    @unpack coeffs = ps
    @views y .= coeffs[:, end]
    @inbounds for i in maxdeg(ps):-1:1
        @views coeffsi = coeffs[:, i]
        @. y = y * x + coeffsi
    end
    y
end

"""
    evaluate!(y::AbstractMatrix, ps::PolySet, x::AbstractVector)

Evaluates a set of polynomials represented by `ps` at multiple scalar points given in `x`, writing the result in-place into the matrix `y`.

# Arguments
- `y`: Output matrix (modified in-place), of size `(n, m)` where `n` is the number of polynomials and `m == length(x)`. Column `j` of `y` contains the evaluations of all polynomials at `x[j]`.
- `ps`: A `PolySet` object storing polynomial coefficients. Each row is a polynomial; columns represent increasing powers of `x`.
- `x`: A vector of scalar inputs at which the polynomials are evaluated.

# Notes
- Uses Horner's method for fast evaluation.
- Evaluation is vectorized over all input points.
- No memory allocation occurs during evaluation.
"""
function evaluate!(y::AbstractMatrix, ps::PolySet, x::AbstractVector)
    @unpack coeffs = ps
    fill!(y, 0)
    @views pend = coeffs[:, end]
    y .+= pend
    @inbounds for i in maxdeg(ps):-1:1
        @views coeffsi = coeffs[:, i]
        @. y = y * x' + coeffsi
    end
    y
end

"""
    integrate!(y::AbstractVector, ps::PolySet{TPS}, a::Real, b::Real) where TPS

Evaluate the definite integral of each polynomial in the `PolySet` `ps` over the interval [`a`, `b`],
and store the result in the output vector `y`.

The result `y[i]` is the integral of the `i`-th polynomial from `a` to `b`.

This version allocates a temporary vector internally.

# Arguments
- `y`: Output vector, must have length equal to the number of polynomials in `ps`.
- `ps`: A `PolySet` structure containing the coefficient matrix.
- `a`, `b`: Integration bounds.
"""
function integrate!(y::AbstractVecOrMat, ps::PolySet{TPS}, a::Real, b::Real) where {TPS}
    @unpack coeffs = ps
    (_, m) = size(ps)
    newT = promote_type(TPS, typeof(a), typeof(b))
    if a != -b
        vec = zeros(newT, m)
        apow = one(newT)
        bpow = one(newT)
        @inbounds for j in 1:m
            apow *= a
            bpow *= b
            vec[j] = (bpow - apow) / j
        end
        mul!(y, coeffs, vec)
        return y
    else
        s = div(m+1, 2)
        vec = zeros(newT, s)
        bpow = one(newT)
        @inbounds for j in 1:s
            bpow *= b
            vec[j] = 2*bpow / (2*j - 1)
            bpow *= b
        end
        @views coeffsv = coeffs[:, 1:2:m]
        mul!(y, coeffsv, vec)
        return y
    end
end

"""
    integrate!(y::AbstractVector, ps::PolySet{TPS}, a::Real, b::Real, cache::AbstractVector) where TPS

Same as `integrate!(y, ps, a, b)` but avoids allocation by reusing the provided `cache` vector.

# Arguments
- `y`: Output vector, must have length equal to the number of polynomials in `ps`.
- `ps`: A `PolySet` structure containing the coefficient matrix.
- `a`, `b`: Integration bounds.
- `cache`: Pre-allocated vector of length equal to the polynomial degree, used to store powers.

"""
function integrate!(y::AbstractVecOrMat, ps::PolySet{TPS}, a::Real,
        b::Real, cache::AbstractVector) where {TPS}
    @unpack coeffs = ps
    (_, m) = size(ps)
    newT = promote_type(TPS, typeof(a), typeof(b))
    apow = one(newT)
    bpow = one(newT)
    @inbounds for j in 1:m
        apow *= a
        bpow *= b
        cache[j] = (bpow - apow) / j
    end
    mul!(y, coeffs, cache)
    y
end

"""
    integrate!(ips::PolySet, ps::PolySet{TPS}, a::Real) where TPS

Computes the indefinite integral of each polynomial in the `PolySet` `ps`, evaluated to be zero at `x = a`,
and stores the result in `ips`.

# Arguments
- `ips::PolySet`: A preallocated `PolySet` where the result will be stored. It must have the same number of
polynomials as `ps` and one more degree.
- `ps::PolySet{TPS}`: A set of polynomials represented as a `PolySet` with scalar type `TPS`.
- `a::Real`: The point at which the antiderivative is set to zero.

# Returns
- `ips::PolySet`: The modified `ips`, containing the integrated polynomials.

# Example
```julia
ps = PolySet([[1.0, 2.0], [0.0, 3.0]])
ips = allocate_polyset(Float64, 2, 3)
integrate!(ips, ps, 0.0)
```
"""
function integrate!(ips::PolySet{TIPS}, ps::PolySet{TPS}, a::Real) where {TIPS,TPS}
    x = 1 ./ (1:(maxdeg(ps) + 1))
    @views vips = ips.coeffs[:, 2:(maxdeg(ps) + 2)]
    mul!(vips, ps.coeffs, Diagonal(x))
    NewT = promote_type(TIPS,TPS, typeof(a))
    y = zeros(NewT, npolys(ps))
    evaluate!(y, ips, a)
    @views vips = ips.coeffs[:, 1]
    vips .-= y
    ips
end

function fill_upper_diagonals!(A::AbstractArray, vals::Base.AbstractVecOrTuple)
    m, n = size(A)
    for d in 1:(length(vals))
        val = vals[d]
        for i in 1:min(m, n - d + 1)
            A[i, i + d - 1] = val
        end
    end
end

function factorize(ps::PolySet{TP}, racine::T) where {TP, T}
    NewT = promote_type(TP, T)
    y = evaluate(ps, racine)
    idx = findall(iszero, y)
    qs = allocate_polyset(NewT, npolys(ps), maxdeg(ps)-1)
    @views vps = ps.coeffs[idx, :]
    @views vqs = qs.coeffs[idx, :]
    M = zeros(NewT, maxdeg(ps)+1, maxdeg(ps)+1)
    fill_upper_diagonals!(M, [-racine, 1])
    @views vc = (vps / M)[:, 1:(end - 1)]
    vqs .= vc
    qs
end

"""
    mul(ps::PolySet{TP}, qs::PolySet{TQ}) where {TP, TQ}

Returns a new `PolySet` containing the product of every polynomial in `ps` with every polynomial in `qs`.

# Arguments
- `ps::PolySet{TP}`: A set of `np` polynomials, each of degree `degp`.
- `qs::PolySet{TQ}`: A set of `nq` polynomials, each of degree `degq`.

# Returns
- `PolySet{NewT}`: A new `PolySet` of `np × nq` polynomials, each of degree at most `degp + degq`,
   where `NewT` is the promoted type of `TP` and `TQ`.

# Details
This function constructs a temporary matrix `M` used to perform the polynomial multiplication via matrix multiplication.
For each polynomial in `qs`, it fills `M` with shifted copies of its coefficients along diagonals, then computes all products with `ps` in one matrix-matrix multiplication.
"""
function mul(ps::PolySet{TP}, qs::PolySet{TQ}) where {TP, TQ}
    (np, degp) = size(ps)
    (nq, degq) = size(qs)
    NewT = promote_type(TP, TQ)
    result = allocate_polyset(NewT, np*nq, degp+degq-2)
    M = zeros(NewT, degp, degp+degq-1)
    @inbounds for i in 1:nq
        @views coeffs = qs.coeffs[i, :]
        fill_upper_diagonals!(M, coeffs)
        @views vresult = result.coeffs[((i - 1) * np + 1):(i * np), :]
        mul!(vresult, ps.coeffs, M)
    end
    result
end

pairwiseproduct(ps::PolySet) = mul(ps, ps)
