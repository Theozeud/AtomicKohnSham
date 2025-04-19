"""
    struct P1IntLegendreGenerator{T} <: AbstractGenerator{T}

Generator for integrated Legendre polynomials and P1 and their derivatives.

# Attributes
- `polynomials::PolySet{T}` : PolySet containing the integrated Legendre polynomials up to order `ordermax`.
- `derivpolynomials::PolySet{T}` : PolySet containing the derivatives of the Legendre polynomials.
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
struct P1IntLegendreGenerator{T, polyType <: PolySet, derivpolyType <: PolySet} <: AbstractGenerator{T}
    polynomials::polyType
    derivpolynomials::derivpolyType
    size::Int
    ordermax::Int
    binf::T
    bsup::T

    function P1IntLegendreGenerator(T::Type = Float64; ordermax::Int = 2, binf::Real = -T(1), bsup::Real = T(1))
        @assert ordermax ≥ 1
        polynomials = allocate_polyset(T, ordermax+1, ordermax)
        derivpolynomials = allocate_polyset(T, ordermax+1, ordermax-1)
        polynomials.coeffs[1,1] = one(T)
        polynomials.coeffs[1,2] = one(T)
        polynomials.coeffs[2,1] = one(T)
        polynomials.coeffs[2,2] = -one(T)
        derivpolynomials.coeffs[1,1] = one(T)
        derivpolynomials.coeffs[2,1] = -one(T)
        leg = Legendre(ordermax-1,T)
        intleg = IntLegendre(ordermax-1, T)
        @views polynomials.coeffs[3:end,:] .= intleg.coeffs[2:end,:]
        @views derivpolynomials.coeffs[3:end,:] .= leg.coeffs[2:end,:]
        new{T, 
            typeof(polynomials), 
            typeof(derivpolynomials)}(polynomials, 
                                      derivpolynomials, 
                                      ordermax + 1, 
                                      ordermax, 
                                      T(binf), 
                                      T(bsup))
    end
end

@inline Base.firstindex(::P1IntLegendreGenerator) = 1
@inline Base.eachindex(p1ilg::P1IntLegendreGenerator) = eachindex(p1ilg.polynomials)
@inline Base.getindex(p1ilg::P1IntLegendreGenerator, n::Int) =  p1ilg.polynomials[n] 
@inline getpolynomial(p1ilg::P1IntLegendreGenerator) = p1ilg.polynomials
@inline getderivpolynomial(p1ilg::P1IntLegendreGenerator) = p1ilg.derivpolynomials
@inline Base.length(p1ilg::P1IntLegendreGenerator) = length(p1ilg.polynomials)

degmax(p1ilg::P1IntLegendreGenerator) = p1ilg.ordermax

function P1IntLegendreGenerator(mesh::Mesh, T::Type = Float64; kwargs...)
    # CREATE GENERATORS
    generators = P1IntLegendreGenerator(T; kwargs...)
    # SIZE OF THE BASIS
    size = (generators.ordermax - 1)* (lastindex(mesh) - 1) +  (lastindex(mesh) - 2)
    # DICTIONNARY TO HAVE ENOUGHT INFORMATIONS ON THE BASIS TO FILL FEM MATRICES EFFICIENTLY
    indices_cells       = Dict{Int,Union{Int,Tuple{Int,Int}}}() #zeros(Int, size, 2)
    indices_generators  = Dict{Int,Union{Int,Tuple{Int,Int}}}() #zeros(Int, size, 2)
    cells_to_indices    = Dict{Int,Vector{Int}}()               #zeros(Int, length(mesh)-1, generators.size)
    cells_to_generators = Dict{Int,Vector{Int}}()
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
    cells_to_generators[i] = collect(2:generators.ordermax+1)
    cells_to_generators[i][1] = 1
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
        cells_to_generators[i] = collect(1:generators.ordermax+1)
        cells_to_generators[i][1] = 2
        cells_to_generators[i][2] = 1
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
    cells_to_generators[i] = collect(2:generators.ordermax+1)

    PolynomialBasis(generators, 
                    mesh, 
                    size, 
                    indices_cells, 
                    indices_generators, 
                    cells_to_indices, 
                    cells_to_generators) 
end