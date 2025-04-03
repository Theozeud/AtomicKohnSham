mutable struct PiecewiseLaurentPolynomial{T,TM} <: AbstractPolynomial{T}
    mesh::Mesh{TM}
    functions::Vector{LaurentPolynomial{T}}
    index::Vector{Int}
    default_value::T
end

@inline Base.zero(pwlp::PiecewiseLaurentPolynomial{T}) where T = zero_piecewiselaurantpolynomial(pwlp.mesh, T)
@inline Base.zero(::Type{PiecewiseLaurentPolynomial{T}}, mesh::Mesh) where T = zero_piecewiselaurantpolynomial(pwlp.mesh, T)
@inline zero_piecewiselaurantpolynomial(mesh::Mesh, T::Type = Float64) = PiecewiseLaurentPolynomial(mesh, LaurentPolynomial{T}[], Int[], T(0))

@inline Base.eachindex(pwlp::PiecewiseLaurentPolynomial) = eachindex(pwlp.mesh)
@inline Base.firstindex(pwlp::PiecewiseLaurentPolynomial) = firstindex(pwlp.mesh)
@inline Base.lastindex(pwlp::PiecewiseLaurentPolynomial) = lastindex(pwlp.mesh)
@inline Base.getindex(pwlp::PiecewiseLaurentPolynomial, i::Int) = i ∈ pwlp.index ? pwlp.functions[findfirst(item->item == i, pwlp.index)] : x->pwlp.default_value

@inline getmesh(pwlp::PiecewiseLaurentPolynomial, i::Int) = pwlp.mesh[i]
@inline getfunction(pwlp::PiecewiseLaurentPolynomial, i::Int) = i ∈ pwlp.index ? pwlp.functions[findfirst(item->item == i, pwlp.index)] : x->pwlp.default_value

(pwlp::PiecewiseLaurentPolynomial)(x) = pwlp[KohnShamResolution.findindex(pwlp.mesh, x)](x) 


function elag!(p::PiecewiseLaurentPolynomial)
    for i ∈ p.index
        i_fun = findfirst(item->item == i, p.index)
        fp = p.functions[i_fun]
        if iszero(fp)
            remove!(p.index, i)
            remove!(p.functions, i_fun)
        end
    end 
    p
end

##################################################################################
#                            Elementary Computations
##################################################################################

function Base.:+(p::PiecewiseLaurentPolynomial{TP}, x::T) where{TP, T}
    NewT = promote_type(TP,T)
    if iszero(x)
        return p
    end
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    for i∈ p.index        
        fp = getfunction(p,i)
        if !ismonomial(fp) || degmin(fp) != 0 || fp[degmin(fp)] != -x
            push!(laurent_poly, x + fp)
            push!(index, i)
        end
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly,index, NewT(x) + NewT(p.default_value))
end

function Base.:+(x, p::PiecewiseLaurentPolynomial)
    p + x
end

function Base.:-(p::PiecewiseLaurentPolynomial{T}) where T
    laurent_poly = LaurentPolynomial{T}[]
    index = Int[]
    for i ∈ p.index
        fp = getfunction(p,i)
        push!(laurent_poly, -fp)
        push!(index, i)
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly,index, -p.default_value)
end

function Base.:-(x, p::PiecewiseLaurentPolynomial)
    x + (-p)
end

function Base.:-(p::PiecewiseLaurentPolynomial, x)
    p + (-x)
end

function Base.:*(p::PiecewiseLaurentPolynomial{TP}, x::T) where{TP, T}
    NewT = promote_type(TP,T)
    if iszero(x)
        return zero(p)
    end
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    for i ∈ p.index
        fp = getfunction(p,i)
        push!(laurent_poly, x * fp)
        push!(index, i)
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly,index, NewT(x) * NewT(p.default_value))
end

function Base.:*(x, p::PiecewiseLaurentPolynomial) 
    p * x
end

function Base.:/(p::PiecewiseLaurentPolynomial, x) 
    p * x^(-1)
end

function Base.:+(p::LaurentPolynomial{TP}, q::PiecewiseLaurentPolynomial{TQ}) where{TP, TQ}
    NewT = promote_type(TP,TQ)
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    for i ∈ eachindex(q)
        if i ∈ q.index
            fq = getfunction(q,i)
            push!(laurent_poly, p + fq)
            push!(index, i)
        else
            push!(laurent_poly, q.default_value + p)
            push!(index, i)
        end
    end
    PiecewiseLaurentPolynomial(q.mesh, laurent_poly, index, NewT(q.default_value))
end

function Base.:+(q::PiecewiseLaurentPolynomial, p::LaurentPolynomial) 
    p + q
end

function Base.:*(p::LaurentPolynomial{TP}, q::PiecewiseLaurentPolynomial{TQ}) where{TP, TQ}
    NewT = promote_type(TP,TQ)
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    if iszero(p)
        return PiecewiseLaurentPolynomial(q.mesh, LaurentPolynomial{NewT}[], Int[], NewT(0))
    end
    for i ∈ eachindex(q)
        if i ∈ q.index
            fq = getfunction(q,i)
            push!(laurent_poly, p * fq)
            push!(index, i)
        elseif q.default_value ≠ 0
            push!(laurent_poly, q.default_value * p)
            push!(index, i)
        end
    end
    PiecewiseLaurentPolynomial(q.mesh, laurent_poly, index, NewT(q.default_value))
end

function Base.:*(q::PiecewiseLaurentPolynomial, p::LaurentPolynomial)
    p * q
end

function Base.:+(p::PiecewiseLaurentPolynomial{TP}, q::PiecewiseLaurentPolynomial{TQ}) where {TP,TQ}
    if p.mesh.points != q.mesh.points
        @error "We can't add for the moment piecewise laurent polynomial on different meshes." 
    end
    NewT = promote_type(TP,TQ)
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    for i ∈ eachindex(p)
        if i ∈ p.index && i ∈ q.index
            fp = getfunction(p,i)
            fq = getfunction(q,i)
            sumfpfq = fp + fq
            elag!(sumfpfq)
            if !iszero(sumfpfq)
                push!(laurent_poly, sumfpfq)
                push!(index,i)
            end
        elseif i ∈ p.index
            fp = getfunction(p,i)
            push!(laurent_poly,  fp + q.default_value)
            push!(index,i)
        elseif i ∈ q.index
            fq = getfunction(q,i)
            push!(laurent_poly , fq + p.default_value)
            push!(index,i)
        end
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly,index, NewT(p.default_value) + NewT(q.default_value))
end

function Base.:*(p::PiecewiseLaurentPolynomial{TP}, q::PiecewiseLaurentPolynomial{TQ}) where{TP, TQ}
    if p.mesh.points != q.mesh.points
        @error "We can't multiply for the moment piecewise laurent polynomial on different meshes." 
    end
    NewT = promote_type(TP,TQ)
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    for i ∈ eachindex(p)
        if i ∈ p.index && i ∈ q.index
            fp = getfunction(p,i)
            fq = getfunction(q,i)
            prodfpfq = fp * fq
            elag!(prodfpfq)
            if !iszero(prodfpfq)
                push!(laurent_poly, prodfpfq)
                push!(index,i)
            end
        elseif q.default_value ≠ 0 && i ∈ p.index
            fp = getfunction(p,i)
            push!(laurent_poly, fp * q.default_value)
            push!(index,i)
        elseif p.default_value ≠ 0 && i ∈ q.index 
            fq = getfunction(q,i)
            push!(laurent_poly, fq * p.default_value)
            push!(index,i)
        end
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly, index, NewT(p.default_value) * NewT(q.default_value))
end

##################################################################################
#                               Composition
##################################################################################
function Base.:∘(p::PiecewiseLaurentPolynomial{TP}, q::LaurentPolynomial{TQ}) where {TP,TQ}
    NewT = promote_type(TP,TQ)
    laurent_poly = LaurentPolynomial{NewT}[]
    index = Int[]
    for i ∈ eachindex(p)
        if i ∈ p.index
            fp = getfunction(p,i)
            push!(laurent_poly, fp∘q)
            push!(index,i)
        end
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly, index, NewT(p.default_value))
end

##################################################################################
#                            Integration & Derivation
##################################################################################

function get_suppport(p::PiecewiseLaurentPolynomial{T}, a::Real, b::Real) where T
    @assert a≤b
    support = Tuple[]
    index_a = KohnShamResolution.findindex(p.mesh, a)
    index_b =  KohnShamResolution.findindex(p.mesh, b)
    if index_a == index_b
        if index_a ∈ p.index
            return [(index_a, a, b)]
        else
            return support
        end
    else
        if index_a ∉ p.index
            index_a += 1
            while index_a ∉ p.index && index_a < index_b
                index_a += 1
            end
            if index_a == index_b
                if index_a ∈ p.index
                    push!(support, (index_a, p.mesh[index_a], b))
                end
                    return support
            else
                push!(support, (index_a, p.mesh[index_a], p.mesh[index_a+1]))
            end
        else
            push!(support,(index_a, a, p.mesh[index_a+1]))
        end
    end
    if index_b ∉ p.index
        index_b -= 1
        while index_b ∉ p.index
            index_b -= 1
        end
        for i in index_a+1:index_b
            push!(support, (i, p.mesh[i], p.mesh[i+1]))
        end
    else
        for i in index_a+1:index_b-1
            push!(support, (i, p.mesh[i], p.mesh[i+1]))
        end
        push!(support,(index_b,p.mesh[index_b],b))
    end
    support
end

function integrate(p::PiecewiseLaurentPolynomial{T}) where T
    laurent_poly = LaurentPolynomial{T}[]
    for i ∈ p.index
        push!(laurent_poly, integrate(p.functions[findfirst(item->item == i, p.index)]))
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly, p.index, T(0))
end

function integrate!(p::PiecewiseLaurentPolynomial)
    for i ∈ p.index
        integrate!(p.functions[findfirst(item->item == i, p.index)])
    end
    p
end

function integrate2(p::PiecewiseLaurentPolynomial{T}, a::Real, b::Real) where T
    @assert a ≤ b
    support = get_suppport(p::PiecewiseLaurentPolynomial{T}, a::Real, b::Real)
    val = T(0)
    for (i,ai,bi) ∈ support
        pi = getfunction(p, i)
        val += integrate(pi, ai, bi)
    end
    val
end

function integrate(p::PiecewiseLaurentPolynomial{T}, a::Real, b::Real) where T
    @assert a ≤ b
    if a == b
        return T(0)
    end
    index_a = KohnShamResolution.findindex(p.mesh, a)
    index_b =  KohnShamResolution.findindex(p.mesh, b)
    val = T(0)
    if index_a + 1> index_b
        if index_a ∈ p.index
            pa = getfunction(p, index_a)
            val += integrate(pa, a, b)
        elseif p.default_value ≠ 0
            val += (b-a) * p.default_value
        end
    else
        if index_a ∈ p.index
            pa = getfunction(p, index_a)
            val += integrate(pa, a, p.mesh[index_a + 1])
        elseif p.default_value ≠ 0
            val += (p.mesh[index_a + 1]-a) * p.default_value
        end
        for i ∈ index_a+1:index_b-1
            if i ∈ p.index
                pᵢ = getfunction(p, i)
                val += integrate(pᵢ, p.mesh[i], p.mesh[i+1])
            elseif p.default_value ≠ 0
                val += (p.mesh[i+1]- p.mesh[i]) * p.default_value
            end
        end
        if index_b ∈ p.index
            pb = getfunction(p, index_b)
            val += integrate(pb, p.mesh[index_b], b)
        elseif p.default_value ≠ 0
            val += (b - p.mesh[index_b]) * p.default_value
        end
    end
    val
end

function deriv(p::PiecewiseLaurentPolynomial{T}) where T
    laurent_poly = LaurentPolynomial{T}[]
    for i ∈ p.index
        push!(laurent_poly, deriv(p.functions[findfirst(item->item == i, p.index)]))
    end
    PiecewiseLaurentPolynomial(p.mesh, laurent_poly, p.index, T(0))
end

function deriv!(p::PiecewiseLaurentPolynomial)
    for i ∈ p.index
        deriv!(p.functions[findfirst(item->item == i, p.index)])
    end
    p
end

scalar_product(p::PiecewiseLaurentPolynomial, q::PiecewiseLaurentPolynomial, a::Real, b::Real) = integrate(p*q,a,b)
scalar_product(p::PiecewiseLaurentPolynomial, q::PiecewiseLaurentPolynomial, m::Mesh) = integrate(p*q, left(m), right(m))

function normL2(p::PiecewiseLaurentPolynomial{T}, m::Mesh) where T
    int = scalar_product(p,p,m)
    if int > 0
        return sqrt(int)
    else
        sqrt(-int)
    end
end
