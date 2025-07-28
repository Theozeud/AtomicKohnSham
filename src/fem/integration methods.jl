abstract type FEMIntegrationMethod end

struct NoSelectedMethod <: FEMIntegrationMethod end

struct ExactIntegration <: FEMIntegrationMethod end

abstract type QuadratureIntegration <: FEMIntegrationMethod end


struct GaussLegendre{T <: Real} <: FEMIntegrationMethod
    npoints::Int
    # FOR FEM MATRIX COMPUTATIONS
    x::Vector{T}
    w::Vector{T}
    Qgenx::Array{T,3}
    # FOR ENERGY COMPUTATIONS
    y::Vector{T}
    wy::Vector{T}
    # CACHE VARIABLES
    shiftx::Vector{T}
    fx::Vector{T}
    fy::Vector{T}
    function GaussLegendre(basis::FEMBasis, npoints::Int = 100)
        @unpack binf, bsup = basis.generators
        # GENERATION OF GAUSS POINTS AND WEIGHTS
        x, w = gausslegendre(npoints)
        # RESCALING ON [0,Rmax]
        Rmax = last(basis.mesh)
        y   = Rmax/2 .*x .+ Rmax/2
        wy  = Rmax/2 .*w
        # RESCALING ON [binf,bsup]
        @. x = (bsup - binf)/2*x + (bsup + binf)/2
        @. w = (bsup - binf)/2*w
        # CACHE VARIABLES FOR FREE-ALLOCATION COMPUTATIONS
        fx = similar(x)
        fy = similar(fx)
        shiftx = similar(fx)
        # EVALUATE POLYNOMIALS
        Qgenx = evaluate(basis.cache.prodMG, x)
        a = Int(sqrt(size(Qgenx,1)))
        Qgenxreshape = reshape(Qgenx, a, a, size(Qgenx,2))
        new{eltype(x)}(npoints, x, w, Qgenxreshape, y, wy, shiftx, fx, fy)
    end
end
