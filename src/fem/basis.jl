#####################################################################
#                           POLYNOMIAL BASIS
#####################################################################

struct PolynomialBasis{ T<:Real, 
                        generatorsType <: AbstractGenerator, 
                        meshType <: Mesh, 
                        dictType <: Dict} <: Basis
    generators::generatorsType                          # Set of Polynomials used to generate the basis
    mesh::meshType                                      # Mesh
    size::Int                                           # Size of the basis
    indices_cells::dictType                             # Dict index basis -> indices support
    indices_generators::dictType                        # Dict index basis -> indices generators
    cells_to_indices::Dict{Int,Vector{Int}}             # Dict index cells -> indices basis
    normalisation::Vector{T}                            # Coefficients of normalisation
    shifts::Vector{Tuple{T,T}}                          # Translation to each cells of the mesh
    invshifts::Vector{Tuple{T,T}}                       # Inverse of the translation 
    matrix_fill_indices::Vector{CartesianIndex{2}}      # Vector of indices filled in fem matrices
    tensor_fill_indices::Vector{CartesianIndex{3}}      # Vector of indices filled in fem tensors
    max_length_intersection::Tuple{Int,Int}             # Maximal number of common mesh two polynomials of the basis can share

    function PolynomialBasis(generators::AbstractGenerator, 
                             mesh::Mesh, 
                             size::Int, 
                             indices_cells::Dict, 
                             indices_generators::Dict,
                             cells_to_indices::Dict{Int,Vector{Int}}, 
                             normalisation::Vector{T}, 
                             shifts::Vector{Tuple{T,T}}, 
                             invshifts::Vector{Tuple{T,T}}, 
                             matrix_fill_indices::Vector{CartesianIndex{2}}, 
                             tensor_fill_indices::Vector{CartesianIndex{3}},
                             max_length_intersection::Tuple{Int,Int}) where T <: Real
        new{eltype(generators), 
            typeof(generators), 
            typeof(mesh),
            typeof(indices_cells)}( generators, 
                                    mesh, 
                                    size, 
                                    indices_cells, 
                                    indices_generators, 
                                    cells_to_indices, 
                                    normalisation, 
                                    shifts, 
                                    invshifts, 
                                    matrix_fill_indices, 
                                    tensor_fill_indices,
                                    max_length_intersection)
    end

    function PolynomialBasis(generators::AbstractGenerator, 
                             mesh::Mesh, 
                             size::Int, 
                             indices_cells::Dict, 
                             indices_generators::Dict, 
                             cells_to_indices::Dict{Int,Vector{Int}}, 
                             normalisation::Vector{<:Real},
                             max_length_intersection::Tuple{Int,Int})
        
        T = eltype(generators)
        shifts    = Vector{Tuple{T,T}}(undef, length(mesh)-1)
        invshifts = Vector{Tuple{T,T}}(undef, length(mesh)-1)
        for i ∈ eachindex(mesh)[1:end-1]
            shifts[i]    = shift(T, mesh[i], mesh[i+1], generators.binf, generators.bsup)
            invshifts[i] = shift(T, generators.binf, generators.bsup, mesh[i], mesh[i+1])
        end
        
        matrix_fill_indices = CartesianIndex{2}[]
        tensor_fill_indices = CartesianIndex{3}[]
        I = fill(0,max_length_intersection[2])
        @inbounds for i in 1:size
            @inbounds for j in i:size
                if !isdisjoint(indices_cells[i], indices_cells[j])
                    push!(matrix_fill_indices, CartesianIndex(i, j))
                    count = find_intersection!(I, indices_cells[i], indices_cells[j])
                    @views Iv = I[1:count]
                    @inbounds for k ∈ j:size
                        if !isdisjoint(Iv, indices_cells[k])
                            push!(tensor_fill_indices, CartesianIndex(i, j, k))
                        end
                    end
                else
                    break
                end
            end
        end 

        new{eltype(generators), 
            typeof(generators),
            typeof(mesh),
            typeof(indices_cells)}( generators, 
                                    mesh, 
                                    size, 
                                    indices_cells, 
                                    indices_generators, 
                                    cells_to_indices, 
                                    normalisation, 
                                    shifts, 
                                    invshifts, 
                                    matrix_fill_indices, 
                                    tensor_fill_indices,
                                    max_length_intersection)
    end
end

@inline Base.eltype(::PolynomialBasis{T, TB, TM}) where {T,TB,TM} = T


@inline Base.length(pb::PolynomialBasis) = pb.size
@inline Base.eachindex(pb::PolynomialBasis) = 1:pb.size

@inline getgenerator(pb::PolynomialBasis, i::Int, j::Int) = getpolynomial(pb.generators, pb.indices_generators[i][j])
@inline getderivgenerator(pb::PolynomialBasis, i::Int, j::Int)= getderivpolynomial(pb.generators, pb.indices_generators[i][j])
@inline getnormalization(pb::PolynomialBasis, i::Int) = pb.normalisation[i]
@inline getshift(pb::PolynomialBasis, i::Int, j::Int) = pb.shifts[pb.indices_cells[i][j]]
@inline getinvshift(pb::PolynomialBasis, i::Int, j::Int) = pb.invshifts[pb.indices_cells[i][j]]
@inline getmesh(pb::PolynomialBasis,i::Int, j::Int) = (pb.mesh[pb.indices_cells[i][j]], pb.mesh[pb.indices_cells[i][j]]+1)

@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    # Linear function that maps [a,b] to [mᵢ, mᵢ₊₁]
    c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
    c0 = -T(a) * T(c1) + T(mᵢ)
    (c1,c0)
end


struct BasisPrecomputations{T <: Real}
    flags::Svector{3,Bool}
    monom_overdeg1::Matrix{T}
    monom_overdeg2::Matrix{T}
    productgenerators::Dict{Tuple{Int,Int},LaurentPolynomial{T}}

    function BasisPrecomputations(mesh::Mesh{TM}, generators::AbstractGenerator{TG}) where{TM,TG}
        T = promote_type(TM,TG)

        flags = @
        degmax(generators)

        productgenerators = Dict{Tuple{Int,Int},LaurentPolynomial{T}}()

        new{T}()
    end
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
    y *= getnormalization(pb, i)
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
    y *= getnormalization(pb, i)
    return y
end