"""
    _integration_monome_over_deg1(k::Int, A::Real, B::Real, a::Real, b::Real)

Calculates the integral :  ∫ₐᵇ Xᵏ/(AX +B) dX

# Arguments:
- `k::Int` : The exponent of the term `X^k` in the numerator. It must be a non-negative integer (`k ≥ 0`).
- `A::Real` : The coefficient multiplying the term `X` in the denominator (`A ≠ 0` if `B ≠ 0`).
- `B::Real` : The constant coefficient in the denominator (`B ≠ 0` if `A ≠ 0`).
- `a::Real` : The lower limit of the integral.
- `b::Real` : The upper limit of the integral.

The integral calculation is performed differently depending on the values of `A` and `B`:
1. **If `B = 0` or `A = 0`**: Explicit formulas are used.
2. **If `A ≠ 0` and `B ≠ 0`**: The function uses the hypergeometric function `pFq` to evaluate the integral.
"""
function _integration_monome_over_deg1(k::Int, A::Real, B::Real, a::Real, b::Real)
    @assert !iszero(A) || !iszero(B)
    @assert k ≥ 0
    @assert !(a ≤ -B/A ≤ b) "Can not integrate X^($k)/($A X + $B) over ($a,$b)"
    if iszero(B)
        if k == 0
            return 1/A * log(abs(b/a))
        else
            return 1/(A * k) * (b^k - a^k)
        end
    elseif iszero(A)
        return 1/(B * (k+1)) * (b^(k+1) - a^(k+1))
    else
        C = B/A
        right = b^(k+1)*pFq((1, k+1), (k+2,), -b/C)
        left = a^(k+1)*pFq((1, k+1), (k+2,), -a/C)
        return (right-left)/(A*(C*k+C))
    end
end

"""
    _integration_monome_over_deg2(k::Int, A::Real, B::Real, a::Real, b::Real)

Calculates the integral :  ∫ₐᵇ Xᵏ/(AX +B)² dX

# Arguments:
- `k::Int` : The exponent of the term `X^k` in the numerator. It must be a non-negative integer (`k ≥ 0`).
- `A::Real` : The coefficient multiplying the term `X` in the denominator (`A ≠ 0` if `B ≠ 0`).
- `B::Real` : The constant coefficient in the denominator (`B ≠ 0` if `A ≠ 0`).
- `a::Real` : The lower limit of the integral.
- `b::Real` : The upper limit of the integral.

The integral calculation is performed differently depending on the values of `A` and `B`:
1. **If `B = 0` or `A = 0`**: Explicit formulas are used.
2. **If `A ≠ 0` and `B ≠ 0`**: The function uses the hypergeometric function `pFq` to evaluate the integral.
"""
function _integration_monome_over_deg2(k::Int, A::Real, B::Real, a::Real, b::Real)
    @assert !iszero(A) || !iszero(B)
    @assert k ≥ 0
    @assert !(a ≤ -B/A ≤ b)
    if iszero(A)
        return 1/(B^2 * (k+1)) * (b^(k+1) - a^(k+1))
    elseif iszero(B)
        if k == 0
            return 1/A^2 * (1/a - 1/b)
        elseif k == 1
            return 1/A^2 * log(b/a)
        elseif k ≥ 2
            return 1/(A^2 * (k-1)) * (b^(k-1) - a^(k-1))
        end
    else
        rac = -B/A
        right = b^(k+1)*pFq((2, k+1), (k+2,), b/rac)
        left = a^(k+1)*pFq((2, k+1), (k+2,), a/rac)
        return 1/(B^2 * (k+1)) * (right - left)
    end
end
