##########################       PolynomialBasis      ##########################

struct PolynomialBasis{T<:Real, generatorsType <: AbstractGenerator, meshType <: Mesh} <: Basis
    generators::generatorsType                      # Set of Polynomials used to generate the basis
    mesh::meshType                                  # Mesh
    size::Int                                       # Size of the basis
    indices_cells::Matrix{Int}                      # Matrix index basis -> indices support
    indices_generators::Matrix{Int}                 # Matrix index basis -> indices generators
    cells_to_indices::Matrix{Int}                   # Matrix index cells -> indices basis
    normalisation::Vector{T}                        # Coefficients of normalisation
    shifts::Vector{Tuple{T,T}}                      # Translation to each cells of the mesh
    invshifts::Vector{Tuple{T,T}}                   # Inverse of the translation 
    matrix_fill_indices::Vector{CartesianIndex{2}}  # Vector of indices filled in fem matrices
    tensor_fill_indices::Vector{CartesianIndex{3}}  # Vector of indices filled in fem tensors

    function PolynomialBasis(generators::AbstractGenerator, 
                             mesh::Mesh, 
                             size::Int, 
                             indices_cells, 
                             indices_generators,
                             cells_to_indices, 
                             normalisation::Vector{T}, 
                             shifts::Vector{Tuple{T,T}}, 
                             invshifts::Vector{Tuple{T,T}}, 
                             matrix_fill_indices::Vector{CartesianIndex{2}}, 
                             tensor_fill_indices::Vector{CartesianIndex{3}}) where T <: Real
        new{eltype(generators), 
            typeof(generators), 
            typeof(mesh)}(  generators, 
                            mesh, 
                            size, 
                            indices_cells, 
                            indices_generators, 
                            cells_to_indices, 
                            normalisation, 
                            shifts, 
                            invshifts, 
                            matrix_fill_indices, 
                            tensor_fill_indices)
    end

    function PolynomialBasis(generators::AbstractGenerator, 
                             mesh::Mesh, 
                             size::Int, 
                             indices_cells, 
                             indices_generators, 
                             cells_to_indices, 
                             normalisation::Vector{T}) 
        
        T = eltype(generators)
        shifts    = Vector{Tuple{T,T}}(undef, length(mesh)-1)
        invshifts = Vector{Tuple{T,T}}(undef, length(mesh)-1)
        for i ∈ eachindex(mesh)[1:end-1]
            shifts[i]    = shift(T, mesh[i], mesh[i+1], generators.binf, generators.bsup)
            invshifts[i] = shift(T, generators.binf, generators.bsup, mesh[i], mesh[i+1])
        end
        
        sorted_indices = sortperm(indices_cells[:, 1])
        indices_cells .= indices_cells[sorted_indices,:]
        indices_generators .= indices_generators[sorted_indices,:]
        normalisation .= normalisation[sorted_indices]

        matrix_fill_indices = CartesianIndex{2}[]
        @inbounds for i in 1:size
            @inbounds for j in i:size
                S = intersect(indices_cells[i,:], indices_cells[j,:])
                if !isempty(S) && !(S==[0])
                    push!(matrix_fill_indices, CartesianIndex(i, j))
                else
                    break
                end
            end
        end 

        tensor_fill_indices = CartesianIndex{3}[]
        @inbounds for I in matrix_fill_indices
            i = I[1]
            j = I[2]
            S = intersect(indices_cells[i,:], indices_cells[j,:])
            @inbounds for k ∈ j:size
                S2 = intersect(S, indices_cells[k,:])
                if !isempty(S2) && !(S2==[0])
                    push!(tensor_fill_indices, CartesianIndex(i, j, k))
                end
            end
        end 

        new{eltype(generators), 
            typeof(generators)}(generators, 
                                mesh, 
                                size, 
                                indices_cells, 
                                indices_generators, 
                                cells_to_indices, 
                                normalisation, 
                                shifts, 
                                invshifts, 
                                matrix_fill_indices, 
                                tensor_fill_indices)
    end
end

@inline Base.eltype(::PolynomialBasis{T, TB, TM}) where {T,TB,TM} = T


@inline Base.length(pb::PolynomialBasis) = pb.size
@inline Base.eachindex(pb::PolynomialBasis) = 1:pb.size

@inline getgenerator(pb::PolynomialBasis, i::Int, j::Int) = getpolynomial(pb.generators, pb.indices_generators[i,j])
@inline getderivgenerator(pb::PolynomialBasis, i::Int, j::Int)= getderivpolynomial(pb.generators, pb.indices_generators[i,j])
@inline getnormalization(pb::PolynomialBasis, i::Int) = pb.normalisation[i]
@inline getshift(pb::PolynomialBasis, i::Int, j::Int) = pb.shifts[pb.indices_cells[i,j]]
@inline getinvshift(pb::PolynomialBasis, i::Int, j::Int) = pb.invshifts[pb.indices_cells[i,j]]


@inline function shift(T::Type, a::Real, b::Real, mᵢ::Real, mᵢ₊₁::Real)
    # Linear function that maps [a,b] to [mᵢ, mᵢ₊₁]
    c1 = (T(mᵢ₊₁) - T(mᵢ))/(T(b) - T(a))
    c0 = -T(a) * T(c1) + T(mᵢ)
    (c1,c0)
end

########################       Evaluation tools       ########################

function (pb::PolynomialBasis)(i::Int, x)
    localisation_x = findindex(pb.mesh, x)
    newT = promote_type(eltype(pb), typeof(x))
    y = zero(newT)
    if localisation_x ∈ pb.indices_cells[i,:]
        for j ∈ axes(pb.indices_generators,2)
            if !iszero(j) && pb.indices_cells[i,j] == localisation_x
                P = getgenerator(pb, i, j)
                ϕ = getshift(pb, i, j)
                ϕx = ϕ[1]*x + ϕ[2]
                y += P(ϕx)
            end
        end
        y *= getnormalization(pb, i)
    end
    return y
end

function (pb::PolynomialBasis)(coeffs::AbstractVector, x)
    @assert length(coeffs) == pb.size
    T = eltype(pb)
    y = zero(T)
    for i ∈ eachindex(pb)
        y += coeffs[i] * pb(i, x)
    end
    y
end

function eval_derivative(pb::PolynomialBasis, i::Int, x)
    localisation_x = findindex(pb.mesh, x)
    newT = promote_type(eltype(pb), typeof(x))
    y = zero(newT)
    if localisation_x ∈ pb.indices_cells[i,:]
        for j ∈ axes(pb.indices_generators,2)
            if !iszero(j) && pb.indices_cells[i,j] == localisation_x
                P = getderivgenerator(pb, i, j)
                ϕx = ϕ[1]*x + ϕ[2]
                y += P(ϕx)
            end
        end
        y *= getnormalization(pb, i)
    end
    return y
end