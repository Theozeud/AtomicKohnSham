#####################################################################
#                      INTEGRATION METHOD STRUCTURE
#####################################################################
abstract type IntegrationMethod end

struct ExactIntegration <: IntegrationMethod end

struct QuadratureIntegration <: IntegrationMethod end

#####################################################################
#                         WEIGHT STRUCTURE
#####################################################################

abstract type AbstractWeight end

struct NoWeight <: AbstractWeight end   # w(x) = 1     
struct InvX <: AbstractWeight end       # w(x) = 1/x
struct InvX2 <: AbstractWeight end      # w(x) = 1/x^2  

has_singularity(::NoWeight, domain::Tuple{T,T}) where T <: Real = false
has_singularity(::InvX,  domain::Tuple{T,T}) where T <: Real    = domain[1] ≤ 0 ≤ domain[2]
has_singularity(::InvX2, domain::Tuple{T,T}) where T <: Real    = domain[1] ≤ 0 ≤ domain[2]


#####################################################################
#                         INTEGRATION DATA
#####################################################################

"""
    IntegrationData{weightType, T, polynomialType, methodType}

Structure holding all the data necessary to compute the integral

    ∫ₐᵇ w(x) ⋅  ∏ᵢ Pᵢ(ϕ(x)) dx

#where:
- `Pᵢ` are Laurent polynomials,
- `ϕ` is an affine map from `[a, b]` to `[binf, bsup]`,
- `w` is a weight function.

# Fields
- `weight::weightType`: The weight function `w(x)`.
- `P::Tuple{Vararg{LaurentPolynomial{T}}}`: Polynomials `Pᵢ`.
- `ϕ::Tuple{T, T}`: Coefficients `(α, β)` defining the affine map `ϕ(x) = α * x + β`.
- `invϕ::Tuple{T, T}`: Coefficients of the inverse affine map of `ϕ`.
- `a::T`: Lower bound of the integration domain.
- `b::T`: Upper bound of the integration domain.
- `binf::T`: Lower bound of the target interval for `ϕ(x)`.
- `bsup::T`: Upper bound of the target interval for `ϕ(x)`.
- `method::IntegrationMethod `: Method used to perform the computations of the integral.

# Type Parameters
- `weightType`: Type of the weight function `w(x)` (subtype of `AbstractWeight`).
- `T`: Numeric type (typically `Float64` or `BigFloat`).
- `methodType` : Type of the method to perform the computation of the integral.

"""
struct IntegrationData{weightType <: AbstractWeight, T<:Real, polynomialType <: Tuple{Vararg{LaurentPolynomial{<:Real}}}, methodType<:IntegrationMethod}
    weight::weightType              
    P::polynomialType
    ϕ::Tuple{T,T}
    invϕ::Tuple{T,T}
    a::T
    b::T
    binf::T
    bsup::T
    method::methodType

    function IntegrationData(weight::AbstractWeight, 
                             P::Tuple{Vararg{LaurentPolynomial{<:Real}}},
                             ϕ::Tuple{<:Real,<:Real},
                             invϕ::Tuple{<:Real,<:Real},
                             a::Real,
                             b::Real,
                             binf::Real,
                             bsup::Real,
                             method::IntegrationMethod)

        new{typeof(weight),
            typeof(a),
            typeof(P),
            typeof(method)}(weight,
                            P,
                            ϕ,
                            invϕ,
                            a,
                            b,
                            binf,
                            bsup,
                            method)
    end
end


function IntegrationData(intdata::IntegrationData; 
                         weight::AbstractWeight = intdata.weight,
                         P::Tuple{Vararg{LaurentPolynomial{<:Real}}} = intdata.P,
                         ϕ::Tuple{<:Real,<:Real} = intdata.ϕ,
                         invϕ::Tuple{<:Real,<:Real} = intdata.invϕ,
                         a::Real = intdata.a,
                         b::Real = intdata.b,
                         binf::Real = intdata.binf,
                         bsup::Real = intdata.bsup,
                         method::IntegrationMethod = intdata.method)

    IntegrationData(weight,
                    P,
                    ϕ,
                    invϕ,
                    a,
                    b,
                    binf,
                    bsup,
                    method)
end

#####################################################################
#                  SHIFTED WEIGHTED SCALAR PRODUCT
#####################################################################

function swsp(intdata::IntegrationData)

    @unpack weight, a, b, method = intdata

    !has_singularity(weight, (a, b)) || return singularity_swsp(weight, intdata)
    return swsp(method, weight, intdata)
end


#####################################################################
#                         HANDLE SINGULARITY
#####################################################################
function find_polynomial_factor(Q::LaurentPolynomial, A::Real, B::Real)
    !iszero(Q) || return Q 
    @assert !iszero(A) || !iszero(B) "Can not find P such that Q = 0 x P"
    if iszero(A)
        return Q / B
    elseif iszero(B)
        return LaurentPolynomial(Q.coeffs ./ A, degmin(Q)-1)
    else
        M = zeros(eltype(Q), length(Q), length(Q)-1)
        for i ∈ 1:length(Q)-1
            M[i,i] = B
            M[i+1,i] = A
        end
        coeffsP = M \ Q.coeffs
        return LaurentPolynomial(coeffsP, degmin(Q))
    end
end


singularity_swsp(::NoWeight, intdata::IntegrationData) = swsp(intdata)


function singularity_swsp(::InvX, intdata::IntegrationData)

    @unpack P, ϕ, method = intdata

    I = findfirst(p -> iszero(p(ϕ[2])), P)

    @assert !isnothing(I) "Singularity in the integrals"

    R = find_polynomial_factor(P[I], 1, -ϕ[2])
    newP = Base.setindex(P, R, I)
    newintdata = IntegrationData(intdata; weight = NoWeight(), P = newP)
    ϕ[1]*swsp(method, NoWeight(), newintdata)
end


function singularity_swsp(::InvX2, intdata::IntegrationData)

    @unpack P, ϕ, method = intdata

    @assert count(p -> iszero(p(ϕ[2])), P) ≥ 2 "Singularity in the integrals"

    I, K = findfirsttwo(p -> iszero(p(ϕ[2])),P)
    RI = find_polynomial_factor(P[I], 1, -ϕ[2])
    RK = find_polynomial_factor(P[K], 1, -ϕ[2])
    newP = Base.setindex(Base.setindex(P, RI, I), RK, K)
    newintdata = IntegrationData(intdata; weight = NoWeight(), P = newP)
    ϕ[1]^2*swsp(method, NoWeight(), newintdata)
end

#####################################################################
#               EXACT SHIFTED WEIGHTED SCALAR PRODUCT
#####################################################################

@memoize memoize_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, a::Real, b::Real) = scalar_product(p, q, a, b)
@memoize memoize_scalar_product(p::LaurentPolynomial, q::LaurentPolynomial, l::LaurentPolynomial, a::Real, b::Real) = scalar_product(p, q*l, a, b)

function swsp(::ExactIntegration, ::NoWeight, intdata::IntegrationData)
    @unpack P, binf, bsup, invϕ = intdata
    invϕ[1] * memoize_scalar_product(P..., binf, bsup)
end

function swsp(::ExactIntegration, ::InvX, intdata::IntegrationData)
    @unpack P, invϕ, binf, bsup = intdata
    Q = reduce(*, P)
    val = zero(eltype(Q))
    for k ∈ eachindex(Q)
        val += Q[k] * _integration_monome_over_deg1(k, invϕ[1], invϕ[2], binf, bsup)
    end 
    invϕ[1] * val
end

function swsp(::ExactIntegration, ::InvX2, intdata::IntegrationData)
    @unpack P, invϕ, binf, bsup = intdata
    Q = reduce(*, P)
    val = zero(eltype(Q))
    for k ∈ eachindex(Q)
        val += Q[k] * _integration_monome_over_deg2(k, invϕ[1], invϕ[2], binf, bsup)
    end 
    invϕ[1] * val
end