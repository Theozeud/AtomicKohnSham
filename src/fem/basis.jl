#####################################################################
#                               PRECOMPUTATIONS
#####################################################################

struct BasisCache{  prodMType <: PolySet,
                    prodTType <: PolySet,
                    T <: Real}
    prodMG::prodMType
    prodMdG::prodMType
    prodTG::prodTType
    K::Matrix{T}                    # Local FEM Matrix
    T::Array{T,3}                   # Local FEM Tensor
    cacheM::Vector{T}
    cacheT::Vector{T}
    evalM::Matrix{T}

    function BasisCache(generators::AbstractGenerator{TG}; 
                        tensor_alloc::Bool = false) where {TG}

        @unpack polynomials, derivpolynomials = generators
        
        # ALL PRODUCT PiPj
        prodMG = pairwiseproduct(polynomials)
        # ALL PRODUCT Pi'Pj'
        prodMdG = pairwiseproduct(derivpolynomials)
        # ALL PRODUCT PiPjPk
        prodTG = tensor_alloc ? mul(prodMG, polynomials) : allocate_polyset(TG,0,-1)

        # LOCAL MATRIX/TENSOR
        K = zeros(TG, size(polynomials,1), size(polynomials,1))
        T = zeros(TG, size(polynomials,1), size(polynomials,1), size(polynomials,1))

        # CACHE
        cacheM = zeros(TG, size(prodMG,2))
        cacheT = zeros(TG, size(prodTG,2))
        evalM = zeros(TG,1,1)

        new{typeof(prodMG),
            typeof(prodTG),
            TG}(prodMG, 
                prodMdG, 
                prodTG,
                K,
                T,
                cacheM,
                cacheT,
                evalM)
    end
end


function getcache(cache::BasisCache, s::Symbol)
    if s == :M
        return cache.cacheM
    elseif s == :Md
        return @views cache.cacheM[1:end-1]
    elseif s == :T
        return cache.cacheT
    end
end


"""
    struct PolynomialBasis{T<:Real, 
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
struct PolynomialBasis{ T<:Real, 
                        generatorsType <: AbstractGenerator, 
                        meshType <: Mesh, 
                        dictType <: Dict,
                        cacheType <: BasisCache} <: Basis
    generators::generatorsType                         
    mesh::meshType                                      
    size::Int                                          
    indices_cells::dictType                            
    indices_generators::dictType                        
    cells_to_indices::Dict{Int,Vector{Int}}             
    cells_to_generators::Dict{Int,Vector{Int}}           
    shifts::Vector{Tuple{T,T}}                          
    invshifts::Vector{Tuple{T,T}}                       
    cache::cacheType                                    

    function PolynomialBasis(generators::AbstractGenerator, 
                             mesh::Mesh, 
                             size::Int, 
                             indices_cells::Dict, 
                             indices_generators::Dict, 
                             cells_to_indices::Dict{Int,Vector{Int}},
                             cells_to_generators::Dict{Int,Vector{Int}})
        
        T = eltype(generators)
        shifts    = Vector{Tuple{T,T}}(undef, length(mesh)-1)
        invshifts = Vector{Tuple{T,T}}(undef, length(mesh)-1)
        for i ∈ eachindex(mesh)[1:end-1]
            shifts[i]    = shift(T, mesh[i], mesh[i+1], generators.binf, generators.bsup)
            invshifts[i] = shift(T, generators.binf, generators.bsup, mesh[i], mesh[i+1])
        end
        
        cache = BasisCache(generators; tensor_alloc = true)

        new{eltype(generators), 
            typeof(generators),
            typeof(mesh),
            typeof(indices_cells),
            typeof(cache)}( generators, 
                            mesh, 
                            size, 
                            indices_cells, 
                            indices_generators, 
                            cells_to_indices, 
                            cells_to_generators,
                            shifts, 
                            invshifts, 
                            cache)
    end
end

@inline Base.eltype(::PolynomialBasis{T}) where {T} = T
@inline Base.length(pb::PolynomialBasis) = pb.size
@inline Base.eachindex(pb::PolynomialBasis) = 1:pb.size

@inline getgenerator(pb::PolynomialBasis, i::Int, j::Int) = getpolynomial(pb.generators, pb.indices_generators[i][j])
@inline getderivgenerator(pb::PolynomialBasis, i::Int, j::Int)= getderivpolynomial(pb.generators, pb.indices_generators[i][j])
@inline getshift(pb::PolynomialBasis, i::Int, j::Int) = pb.shifts[pb.indices_cells[i][j]]
@inline getinvshift(pb::PolynomialBasis, i::Int, j::Int) = pb.invshifts[pb.indices_cells[i][j]]
@inline getmesh(pb::PolynomialBasis,i::Int, j::Int) = (pb.mesh[pb.indices_cells[i][j]], pb.mesh[pb.indices_cells[i][j]]+1)


function Base.show(io::IO, basis::PolynomialBasis)
    println(io, "PolynomialBasis with $(basis.size) functions")
    println(io, "  Mesh type:        $(typeof(basis.mesh))")
    println(io, "  Generators type:  $(typeof(basis.generators))")
    println(io, "  Cells:            $(length(basis.shifts)) total")
    println(io, "  Caching:          $(typeof(basis.cache))")
end


@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    # Linear function that maps [a,b] to [mᵢ, mᵢ₊₁]
    c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
    c0 = -T(a) * T(c1) + T(mᵢ)
    (c1,c0)
end


function getelement(pb::PolynomialBasis, k::Int, s::Symbol)
    @unpack generators, mesh, shifts, invshifts = pb
    ElementData(shifts[k], 
                invshifts[k], 
                mesh[k], 
                mesh[k+1], 
                generators.binf, 
                generators.bsup,
                s)
end

#####################################################################
#                          EVALUATION TOOLS
#####################################################################
function (pb::PolynomialBasis)(i::Int, x::T) where T
    localisation_x = findindex(pb.mesh, x)
    j = findfirst(==(localisation_x), pb.indices_cells[i])
    if isnothing(j)
        newT = promote_type(eltype(pb), T)
        return zero(newT)
    end
    P = getgenerator(pb, i, j)
    ϕ = getshift(pb, i, j)
    ϕx = ϕ[1]*x + ϕ[2]
    y = P(ϕx)
    return y
end

function (pb::PolynomialBasis)(coeffs::AbstractVector, x::T) where T
    @assert length(coeffs) == pb.size
    newT = promote_type(eltype(pb),T)
    y = zero(newT)
    localisation_x = findindex(pb.mesh, x)
    for i ∈ cells_to_indices[localisation_x]
        y += coeffs[i] * pb(i, x)
    end
    y
end

function eval_derivative(pb::PolynomialBasis, i::Int, x::T) where T
    localisation_x = findindex(pb.mesh, x)
    j = findfirst(==(localisation_x), pb.indices_cells[i])
    if isnothing(j)
        newT = promote_type(eltype(pb), T)
        return zero(newT)
    end
    P = getderivgenerator(pb, i, j)
    ϕ = getshift(pb, i, j)
    ϕx = ϕ[1]*x + ϕ[2]
    y = P(ϕx)
    return y
end