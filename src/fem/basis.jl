#--------------------------------------------------------------------
#                               CACHE
#--------------------------------------------------------------------
struct BasisCache{prodMType <: PolySet,
    prodTType <: PolySet,
    T <: Real}
    prodMG::prodMType
    prodMdG::prodMType
    prodTG::prodTType
    K::Matrix{T}                    # Local FEM Matrix
    T::Array{T, 3}                   # Local FEM Tensor
    cacheM::Vector{T}
    cacheT::Vector{T}

    function BasisCache(generators::AbstractGenerator{TG}) where {TG}
        @unpack polynomials, derivpolynomials = generators

        # ALL PRODUCT PiPj
        prodMG = pairwiseproduct(polynomials)
        # ALL PRODUCT Pi'Pj'
        prodMdG = pairwiseproduct(derivpolynomials)
        # ALL PRODUCT PiPjPk
        prodTG = mul(prodMG, polynomials)

        # LOCAL MATRIX/TENSOR
        K = zeros(TG, size(polynomials, 1), size(polynomials, 1))
        T = zeros(TG, size(polynomials, 1), size(polynomials, 1), size(polynomials, 1))

        # CACHE
        cacheM = zeros(TG, size(prodMG, 2))
        cacheT = zeros(TG, size(prodTG, 2))

        new{typeof(prodMG),
            typeof(prodTG),
            TG}(prodMG,
            prodMdG,
            prodTG,
            K,
            T,
            cacheM,
            cacheT)
    end
end

function getcache(cache::BasisCache, s::Symbol)
    if s == :M
        return cache.cacheM
    elseif s == :Md
        return @views cache.cacheM[1:(end - 1)]
    elseif s == :T
        return cache.cacheT
    end
end

#--------------------------------------------------------------------
# 						FEM BASIS STRUCTURE
#--------------------------------------------------------------------
"""
    struct FEMBasis{T<:Real, 
                           generatorsType <: AbstractGenerator, 
                           meshType <: Mesh, 
                           dictType <: Dict,
                           cacheType <: BasisCache} <: Basis

Represents a basis of polynomial functions defined over a 1D mesh. 
The basis is constructed from a set of polynomial generators and adapted locally to each cell of the mesh.

# Fields

- `generators::generatorsType`  
  Set of polynomial generators used to build the basis. These can be Legendre polynomials, monomials, or any other suitable family.

- `mesh::meshType`  
  The mesh over which the basis is defined. 

- `size::Int`  
  Total number of basis functions.

- `indices_cells::dictType`  
  Dictionary mapping each basis index to the cell(s) on which its support is nonzero.

- `indices_generators::dictType`  
  Dictionary mapping each basis index to the generator indices used to construct it.

- `cells_to_indices::Dict{Int,Vector{Int}}`  
  Inverse map: associates each cell with the indices of basis functions that are supported on it.

- `cells_to_generator::Dict{Int,Vector{Int}}`  
  Maps each cell to the generator indices used locally for constructing basis functions on that cell.

- `shifts::Vector{Tuple{T,T}}`  
  Translations applied to center each cell in a reference domain (e.g., mapping to a reference element).

- `invshifts::Vector{Tuple{T,T}}`  
  Inverse translations used to map from the reference cell back to the global coordinates.

- `cache::cacheType`  
  A cache structure used to store intermediate values for efficient computation.

"""
struct FEMBasis{T <: Real,
    generatorsType <: AbstractGenerator,
    meshType <: Mesh,
    dictType <: Dict,
    cacheType <: BasisCache}
    generators::generatorsType
    mesh::meshType
    size::Int
    indices_cells::dictType
    indices_generators::dictType
    cells_to_indices::Dict{Int, Vector{Int}}
    cells_to_generators::Dict{Int, Vector{Int}}
    shifts::Vector{Tuple{T, T}}
    invshifts::Vector{Tuple{T, T}}
    cache::cacheType
    max_nb_poly_cells::Int

    function FEMBasis(generators::AbstractGenerator,
            mesh::Mesh,
            size::Int,
            indices_cells::Dict,
            indices_generators::Dict,
            cells_to_indices::Dict{Int, Vector{Int}},
            cells_to_generators::Dict{Int, Vector{Int}})
        T = eltype(generators)
        shifts = Vector{Tuple{T, T}}(undef, length(mesh)-1)
        invshifts = Vector{Tuple{T, T}}(undef, length(mesh)-1)
        for i in eachindex(mesh)[1:(end - 1)]
            shifts[i] = shift(T, mesh[i], mesh[i + 1], generators.binf, generators.bsup)
            invshifts[i] = shift(T, generators.binf, generators.bsup, mesh[i], mesh[i + 1])
        end

        cache = BasisCache(generators)

        max_nb_poly_cells = max(length.(values(cells_to_indices))...)

        new{eltype(generators),
            typeof(generators),
            typeof(mesh),
            typeof(indices_cells),
            typeof(cache)}(generators,
            mesh,
            size,
            indices_cells,
            indices_generators,
            cells_to_indices,
            cells_to_generators,
            shifts,
            invshifts,
            cache,
            max_nb_poly_cells)
    end
end

#--------------------------------------------------------------------
#								API
#--------------------------------------------------------------------
@inline Base.eltype(::FEMBasis{T}) where {T} = T
@inline Base.length(pb::FEMBasis) = pb.size
@inline Base.eachindex(pb::FEMBasis) = 1:pb.size

@inline getgenerator(pb::FEMBasis, i::Int,
    j::Int) = getpolynomial(pb.generators, pb.indices_generators[i][j])
@inline getderivgenerator(pb::FEMBasis, i::Int,
    j::Int) = getderivpolynomial(pb.generators, pb.indices_generators[i][j])
@inline getshift(pb::FEMBasis, i::Int, j::Int) = pb.shifts[pb.indices_cells[i][j]]
@inline getinvshift(pb::FEMBasis, i::Int, j::Int) = pb.invshifts[pb.indices_cells[i][j]]
@inline getmesh(pb::FEMBasis, i::Int,
    j::Int) = (pb.mesh[pb.indices_cells[i][j]], pb.mesh[pb.indices_cells[i][j]]+1)

function Base.show(io::IO, basis::FEMBasis)
    println(io, "FEMBasis with $(basis.size) functions")
    println(io, "  Mesh type:        $(typeof(basis.mesh))")
    println(io, "  Generators type:  $(typeof(basis.generators))")
    println(io, "  Cells:            $(length(basis.shifts)) total")
    println(io, "  Caching:          $(typeof(basis.cache))")
end

@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    # Linear function that maps [a,b] to [mᵢ, mᵢ₊₁]
    c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
    c0 = -T(a) * T(c1) + T(mᵢ)
    (c1, c0)
end

#--------------------------------------------------------------------
#                          EVALUATION TOOLS
#--------------------------------------------------------------------
```
    Evaluation of the i-th element of the basis `pb` in x.

```
function (pb::FEMBasis)(i::Int, x::T) where {T}
    localisation_x = findindex(pb.mesh, x)
    j = findfirst(==(localisation_x), pb.indices_cells[i])
    if isnothing(j)
        newT = promote_type(eltype(pb), T)
        return zero(newT)
    end
    P = getgenerator(pb, i, j)
    ϕ = getshift(pb, i, j)
    ϕx = ϕ[1]*x + ϕ[2]
    y = evaluate(P, ϕx)
    return y
end

```
    Evaluation of the i-th element of the basis `pb` in x
    for i ∈ `I`.

```
function (pb::FEMBasis)(I::AbstractVector{Int}, x::T) where {T}
    NewT = promote_type(eltype(pb), T)
    evaluations = zeros(NewT, length(I))
    evaluate!(evaluations, pb, I, x)
    evaluations
end

function (pb::FEMBasis)(x::T) where {T}
    NewT = promote_type(eltype(pb), T)
    evaluations = zeros(NewT, length(I))
    evaluate!(evaluations, pb, 1:length(pb), x)
    evaluations
end

function evaluate!(evaluations::AbstractVector,
        pb::FEMBasis,
        I::AbstractVector{<:Int},
        x::T,
        cache_Pϕx::AbstractVector = _cache_Pϕx(pb, x)) where {T}
    fill!(evaluations, 0)
    fill!(cache_Pϕx, 0)
    localisation_x = findindex(pb.mesh, x)
    ϕ = pb.shifts[localisation_x]
    ϕx = ϕ[1]*x + ϕ[2]
    P = pb.generators.polynomials
    AtomicKohnSham.evaluate!(cache_Pϕx, P, ϕx)
    cache_Pϕx
    @inbounds for i in eachindex(I)
        j = findfirst(==(localisation_x), pb.indices_cells[I[i]])
        if !isnothing(j)
            k = pb.indices_generators[I[i]][j]
            evaluations[i] = cache_Pϕx[k]
        end
    end
end

function _cache_Pϕx(pb::FEMBasis, ::T) where {T}
    newT = promote_type(eltype(pb), T)
    zeros(newT, length(pb.generators))
end

```
    Evaluation of a function defined through coefficients `coeffs` over the basis `pb` 
	at the points of the vector `X`. 

```
function evaluate(
        pb::FEMBasis, coeffs::AbstractVector{<:Real}, X::AbstractVector{T}) where {T}
    NewT = promote_type(eltype(pb), eltype(coeffs), T)
    evaluations = zeros(NewT, length(X))
    evaluate!(evaluations, pb, coeffs, X)
    evaluations
end

function evaluate!(evaluations::AbstractVector,
        pb::FEMBasis,
        coeffs::AbstractVector{<:Real},
        X::AbstractVector{<:Real})
    fill!(evaluations, 0)
    Q = zeros(eltype(evaluations), pb.max_nb_poly_cells)
    cache_Pϕx = _cache_Pϕx(pb, first(X))
    @inbounds for k in eachindex(X)
        xk = X[k]
        localisation_xk = findindex(pb.mesh, xk)
        Ik = pb.cells_to_indices[localisation_xk]
        lk = length(Ik)
        @views coeffsk = coeffs[Ik]
        @views Qk = Q[1:lk]
        evaluate!(Qk, pb, Ik, xk, cache_Pϕx)
        evaluations[k] = dot(Qk, coeffsk)
    end
end
