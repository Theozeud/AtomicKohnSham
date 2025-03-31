"""
    struct P1IntLegendreGenerator{T} <: AbstractGenerator{T}

Generator for integrated Legendre polynomials and P1 and their derivatives.

# Attributes
- `polynomials::Vector{LaurentPolynomial{T}}` : Vector containing the integrated Legendre polynomials up to order `ordermax`.
- `derivpolynomials::Vector{LaurentPolynomial{T}}` : Vector containing the derivatives of the Legendre polynomials.
- `size::Int` : Size of the polynomial vector, equal to `ordermax + 1`.
- `ordermax::Int` : Maximum order of the generated polynomials.
- `binf::T` : Lower bound of the definition interval.
- `bsup::T` : Upper bound of the definition interval.

# Constructor
    P1IntLegendreGenerator(T::Type = Float64; ordermax = 2, binf::Real = -T(1), bsup::Real = T(1))

Constructs a generator for integrated Legendre polynomials.

## Arguments
- `T::Type` : Numeric type for polynomial coefficients (default `Float64`).
- `ordermax::Int` : Maximum polynomial order (must be `≥ 1`).
- `binf::Real` : Lower bound of the definition interval (default `-1`).
- `bsup::Real` : Upper bound of the definition interval (default `1`).
"""
struct P1IntLegendreGenerator{T} <: AbstractGenerator{T}
    polynomials::Vector{LaurentPolynomial{T}}
    derivpolynomials::Vector{LaurentPolynomial{T}}
    size::Int
    ordermax::Int
    binf::T
    bsup::T

    function P1IntLegendreGenerator(T::Type = Float64; ordermax::Int = 2, binf::Real = -T(1), bsup::Real = T(1))
        @assert ordermax ≥ 1
        polynomials = Vector{LaurentPolynomial{T}}(undef, ordermax+1)
        derivpolynomials = Vector{LaurentPolynomial{T}}(undef, ordermax+1)
        polynomials[1] = Polynomial([one(T),one(T)], 0)
        polynomials[2] = Polynomial([one(T),-one(T)], 0)
        derivpolynomials[1] = Polynomial([one(T)], 0)
        derivpolynomials[2] = Polynomial([-one(T)], 0)
        for n ∈ 2:ordermax
            Pₙ = Legendre(n-1; T = T, a = T(binf), b = T(bsup))
            derivpolynomials[n+1] = Pₙ
            Qₙ = intLegendre(n-1; T = T, a = T(binf), b = T(bsup))
            polynomials[n+1] = Qₙ
        end
        new{T}(polynomials, derivpolynomials, ordermax + 1, ordermax, T(binf), T(bsup))
    end
end

@inline Base.firstindex(::P1IntLegendreGenerator) = 1
@inline Base.eachindex(p1ilg::P1IntLegendreGenerator) = eachindex(p1ilg.polynomials)
@inline Base.getindex(p1ilg::P1IntLegendreGenerator, n::Int) =  p1ilg.polynomials[n] 
@inline getpolynomial(p1ilg::P1IntLegendreGenerator) = p1ilg.polynomials
@inline getderivpolynomial(p1ilg::P1IntLegendreGenerator) = p1ilg.derivpolynomials


function P1IntLegendreGenerator(mesh::Mesh, T::Type = Float64; kwargs...)
    # CREATE GENERATORS
    generators = P1IntLegendreGenerator(T; kwargs...)
    # SIZE OF THE BASIS
    size = (generators.ordermax - 1)* (lastindex(mesh) - 1) +  (lastindex(mesh) - 2)
    # DICTIONNARY TO HAVE ENOUGHT INFORMATIONS ON THE BASIS TO FILL FEM MATRICES EFFICIENTLY
    indices_cells       = Dict{Int,Union{Int,Tuple{Int,Int}}}() #zeros(Int, size, 2)
    indices_generators  = Dict{Int,Union{Int,Tuple{Int,Int}}}() #zeros(Int, size, 2)
    cells_to_indices    = Dict{Int,Vector{Int}}()               #zeros(Int, length(mesh)-1, generators.size)

    numbas = 1
    # First mesh
    i = firstindex(mesh)
    indices_cells[numbas]         = (i,i+1) 
    indices_generators[numbas]    = (1,2)
    numbas += 1
    for j ∈ 2:generators.ordermax
        indices_cells[numbas]       = i
        indices_generators[numbas]  = j+1
        numbas += 1
    end
    cells_to_indices[i] = 1:generators.ordermax
    # Mid-meshes
    for i ∈ eachindex(mesh)[2:end-2]
        indices_cells[numbas]       = (i,i+1) 
        indices_generators[numbas]  = (1,2)
        numbas += 1
        for j ∈ 2:generators.ordermax
            indices_cells[numbas]       = i
            indices_generators[numbas]  = j+1
            numbas += 1
        end
        cells_to_indices[i] = zeros(Int,generators.ordermax+1)
        cells_to_indices[i][1] = numbas - 2*generators.ordermax
        cells_to_indices[i][2:end] .= (numbas - generators.ordermax):numbas-1
    end
    # Last mesh
    i = lastindex(mesh)-1
    for j ∈ 2:generators.ordermax
        indices_cells[numbas]       = i
        indices_generators[numbas]  = j+1
        numbas += 1
    end
    cells_to_indices[i] = zeros(Int,generators.ordermax)
    cells_to_indices[i][1] = numbas - 2*generators.ordermax + 1
    cells_to_indices[i][2:end] .= (numbas - generators.ordermax+1):numbas-1

    # NORMALISATION COEFFICIENTS
    normalisation = ones(Int, size)

    PolynomialBasis(generators, 
                    mesh, 
                    size, 
                    indices_cells, 
                    indices_generators, 
                    cells_to_indices, 
                    normalisation,
                    (2,2)) 
end