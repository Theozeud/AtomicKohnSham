struct Polynomial{typeData<:AbstractVector}
    coeffs::typeData
end

Base.eltype(p::Polynomial) = eltype(p.coeffs)
Base.getindex(p::Polynomial, i::Int) = p.coeffs[i]

function (p::Polynomial)(x::T) where T
    NewT = promote_type(T, eltype(p))
    y = zero(NewT)
    y = a[end]
    for i ∈ length(a)-1:-1:1
        y = y*x + a[i]
    end
    y
end


@memoize integal_monom(k::Int, a::Real, b::Real) = 1/(k+1) * (b^(k+1) - a^(k+1))

