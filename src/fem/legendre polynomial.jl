"""
    Legendre(n::Int, ::Type{T}) where {T}

Returns a matrix `polys` of size `(n+1, n+1)` containing the coefficients of the Legendre polynomials `P₀(x)` to `Pₙ(x)` expressed in the monomial basis `{1, x, x², ..., xⁿ}`.

Each row `i` corresponds to the coefficients of the polynomial `P_{i-1}(x)`, ordered by increasing degree.

# Arguments
- `n::Int`: The maximum degree of the Legendre polynomial to compute.
- `T::Type`: The numerical type used for the coefficients (e.g., `Float64`, `Double64`, etc.).

# Returns
- `Matrix{T}`: A matrix where `polys[i, j]` is the coefficient of `x^(j-1)` in `P_{i-1}(x)`.

# Details
- The first polynomials (`P₀`, `P₁`, `P₂`) are stored in the constant `LEGENDRE_COEFFS` to avoid recomputation.
- For `n > LEGENDRE_DEGMAX_STORED`, the remaining polynomials are computed recursively using the Bonnet recurrence relation:
  `Pₙ(x) = ((2n−1)/n) * x * Pₙ₋₁(x) − ((n−1)/n) * Pₙ₋₂(x)`

# Examples
```julia
Legendre(2)           # Returns coefficients for P₀, P₁, and P₂ as Float64
Legendre(5, Double64) # Returns coefficients for P₀ to P₅ using Double64 precision
```
"""
Legendre(n::Int) = Legendre(n, Float64)

const LEGENDRE_DEGMAX_STORED = 2

const LEGENDRE_COEFFS = [1//1 0//1 0//1;
                        0//1 1//1 0//1;
                        -1//2 0//1 3//2;]

function Legendre(n::Int, ::Type{T}) where {T}
    if n ≤ LEGENDRE_DEGMAX_STORED
        @views polys = T.(LEGENDRE_COEFFS[1:n+1,1:n+1])
        return PolySet(polys)
    else
        polys = zeros(T,n+1,n+1)
        @views vpolys = polys[1:LEGENDRE_DEGMAX_STORED+1,1:LEGENDRE_DEGMAX_STORED+1]
        vpolys .= LEGENDRE_COEFFS
        for m ∈ LEGENDRE_DEGMAX_STORED+1:n
            @views Pₘ₋₁ = polys[m,1:m]
            @views Pₘ₋₂ = polys[m-1,1:m-1]
            @views polys[m+1,2:m+1]  .= (2*m-1)//m .*  Pₘ₋₁ 
            @views polys[m+1,1:m-1] .-= (m - 1)//m .* Pₘ₋₂
        end
        return PolySet(polys)
    end
end


"""
    IntLegendre(n::Int, ::Type{T}) where {T}

Returns the set of indefinite integrals of Legendre polynomials of degrees 0 to `n`, evaluated so that they vanish at `x = -1`. The result is returned as a `PolySet` with coefficients of type `T`.

# Arguments
- `n::Int`: The maximum degree of Legendre polynomials to integrate.
- `T`: The desired scalar type for the coefficients (e.g., `Float64`, `Double64`, etc.).

# Returns
- `PolySet`: A `PolySet` containing the integrals..

# Details
If `n` is small enough (i.e., `n ≤ INTLEGENDRE_DEGMAX_STORED`), the result is retrieved from a precomputed table of rational coefficients for improved performance and precision.

# Example
```julia
julia> ps = IntLegendre(3, Float64)
PolySet with 4 polynomials of max degree 4
```
"""
IntLegendre(n::Int) = IntLegendre(n, Float64)

const INTLEGENDRE_DEGMAX_STORED = 1

const INTLEGENDRE_COEFFS = [1//1    1//1    0//1;
                            -1//2   0//1    1//2]

function IntLegendre(n::Int, ::Type{T}) where {T}
    if n ≤ INTLEGENDRE_DEGMAX_STORED
        @views polys = T.(INTLEGENDRE_COEFFS[1:n+1,1:n+2])
        return PolySet(polys)
    else
        polys = zeros(T,n+1,n+2)
        @views vpolys = polys[1:INTLEGENDRE_DEGMAX_STORED+1,1:INTLEGENDRE_DEGMAX_STORED+2]
        vpolys .= INTLEGENDRE_COEFFS
        legendre = Legendre(n, Rational{Int64})
        @views vpolys = polys[INTLEGENDRE_DEGMAX_STORED+2:end,:]
        @views vlegendre = legendre.coeffs[INTLEGENDRE_DEGMAX_STORED+2:end,:]
        integrate!(PolySet(vpolys), PolySet(vlegendre), -1)
        return PolySet(polys)
    end
end