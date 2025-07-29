abstract type AbstractWeight end

struct NoWeight <: AbstractWeight end   # w(x) = 1     
struct InvX <: AbstractWeight end       # w(x) = 1/x
struct InvX2 <: AbstractWeight end      # w(x) = 1/x^2  

struct FunWeight{funType <: Base.Function} <: AbstractWeight
    is_vectorized::Bool
    is_inplace::Bool
    f::funType
    function FunWeight(f::F; is_vectorized::Bool, is_inplace::Bool) where {F}
        @assert (is_inplace && is_vectorized) || !is_inplace "If f is in place it must be vectorized too."
        new{typeof(f)}(is_vectorized, is_inplace, f)
    end
end

function (weight::FunWeight)(Y::AbstractVector, X::AbstractVector)
    @unpack is_vectorized, is_inplace, f = weight
    if is_inplace
        f(Y, X)
    elseif is_vectorized
        Y .= f(X)
    else
        @. Y = f(X)
    end
end

default_method(::FEMBasis, ::AbstractWeight) = ExactIntegration()
default_method(::FEMBasis, method::FEMIntegrationMethod, ::AbstractWeight) = method
function default_method(basis::FEMBasis, ::NoSelectedMethod, weight::AbstractWeight)
    default_method(basis, weight)
end
default_method(basis::FEMBasis, ::FunWeight) = GaussLegendre(basis)

has_singularity(::NoWeight, domain::Tuple{T, T}) where {T <: Real} = false
has_singularity(::InvX, domain::Tuple{T, T}) where {T <: Real} = domain[1] ≤ 0 ≤ domain[2]
has_singularity(::InvX2, domain::Tuple{T, T}) where {T <: Real} = domain[1] ≤ 0 ≤ domain[2]
has_singularity(::FunWeight, domain::Tuple{T, T}) where {T <: Real} = false
