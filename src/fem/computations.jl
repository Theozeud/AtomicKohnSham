#####################################################################
#                      INTEGRATION METHOD STRUCTURE
#####################################################################
abstract type IntegrationMethod end

struct NoSelectedMethod <: IntegrationMethod end

struct ExactIntegration <: IntegrationMethod end

struct QuadratureIntegration{T <: Real} <: IntegrationMethod 
    x::Vector{T}
    w::Vector{T}
    shiftx::Vector{T}
    wx::Vector{T}
    function QuadratureIntegration()
        x, w = gausslegendre(100)
        wx = similar(x)
        shiftx = similar(wx)
        new{eltype(x)}(x,w,shiftx,wx)
    end
end



#####################################################################
#                         WEIGHT STRUCTURE
#####################################################################
abstract type AbstractWeight end

struct NoWeight <: AbstractWeight end   # w(x) = 1     
struct InvX <: AbstractWeight end       # w(x) = 1/x
struct InvX2 <: AbstractWeight end      # w(x) = 1/x^2  

struct FunWeight{funType <: Base.Function} <: AbstractWeight 
    f::funType
end

(weight::FunWeight)(args...; kwargs...) = weight.f(args...; kwargs...)


default_method(::AbstractWeight) = ExactIntegration()
default_method(method::IntegrationMethod, ::AbstractWeight) = method
default_method(::NoSelectedMethod, weight::AbstractWeight) = default_method(weight)
default_method(::FunWeight) = QuadratureIntegration()


has_singularity(::NoWeight, domain::Tuple{T,T}) where T <: Real = false
has_singularity(::InvX,  domain::Tuple{T,T}) where T <: Real    = domain[1] ≤ 0 ≤ domain[2]
has_singularity(::InvX2, domain::Tuple{T,T}) where T <: Real    = domain[1] ≤ 0 ≤ domain[2]
has_singularity(::FunWeight, domain::Tuple{T,T}) where T<: Real = false



#####################################################################
#                           ELEMENT DATA
#####################################################################
struct ElementData{ T<:Real}
    ϕ::Tuple{T,T}
    invϕ::Tuple{T,T}
    a::T
    b::T
    binf::T
    bsup::T
    s::Symbol

    function ElementData(ϕ::Tuple{<:Real,<:Real},
                        invϕ::Tuple{<:Real,<:Real},
                        a::Real,
                        b::Real,
                        binf::Real,
                        bsup::Real,
                        s::Symbol)

    new{typeof(a)}( ϕ,
                    invϕ,
                    a,
                    b,
                    binf,
                    bsup,
                    s)
    end
end


function ElementData(eldata::ElementData; 
                     ϕ::Tuple{<:Real,<:Real} = eldata.ϕ,
                     invϕ::Tuple{<:Real,<:Real} = eldata.invϕ,
                     a::Real = eldata.a,
                     b::Real = eldata.b,
                     binf::Real = eldata.binf,
                     bsup::Real = eldata.bsup,
                     s::Symbol = eldata.s)

    ElementData(ϕ,
                invϕ,
                a,
                b,
                binf,
                bsup,
                s)
end