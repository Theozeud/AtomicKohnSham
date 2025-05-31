# Using the Bonnet Recurrence formula to compute the nth, normalized or not, Legendre Polynomial scaled on [a,b]
@memoize function Legendre(n::Int; a::Real = -1, b::Real = 1, T::Type = Float64, normalize::Bool = false)
    if n == 0
        if normalize
            return Polynomial([T(1)]) * sqrt(T(1)/T(b-a))
        else
            return Polynomial([T(1)])
        end
    elseif n == 1
        if normalize
            return Polynomial([T(-b-a)/T(b-a), T(2)/T(b-a)]) * sqrt(T(3)/T(b-a))
        else
            return Polynomial([T(-b-a)/T(b-a), T(2)/T(b-a)])
        end
    else
        coeff = normalize ? sqrt(T(2*n + 1)/T(b-a)) : T(1)
        pₙ₋₁ = Legendre(n-1; a = a, b = b, T = T, normalize = false)
        pₙ₋₂ = Legendre(n-2; a = a, b = b, T = T, normalize = false)
        pₙ = (2*T(n)-1)/T(n) * Polynomial([T(-b-a)/T(b-a), T(2)/T(b-a)]) * pₙ₋₁ - (T(n)-T(1))/T(n) * pₙ₋₂
        return coeff * pₙ 
    end
end

# Integrated Legendre Polynomial
@memoize function intLegendre(n::Int; a::Real = -1, b::Real = 1, T::Type = Float64, normalize::Bool = false)
    pₙ = Legendre(n; a = a, b = b, T = T, normalize = false)
    qₙ = integrate(pₙ)
    qₙ = qₙ - qₙ(a)
    if normalize
        return qₙ / sqrt(scalar_product(qₙ, qₙ, a, b))
    else
        return qₙ
    end
end