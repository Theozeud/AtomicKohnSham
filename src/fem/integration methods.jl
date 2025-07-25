abstract type FEMIntegrationMethod end

struct NoSelectedMethod <: FEMIntegrationMethod end

struct ExactIntegration <: FEMIntegrationMethod end

abstract type QuadratureIntegration <: FEMIntegrationMethod end 


struct GaussLegendre{T <: Real} <: FEMIntegrationMethod 
    npoints::Int
    x::Vector{T}
    w::Vector{T}
    shiftx::Vector{T}
    wx::Vector{T}
    Qgenx::Array{T,3}
    function GaussLegendre(basis::FEMBasis, npoints::Int = 100)
        @unpack binf, bsup = basis.generators
        # GENERATION OF GAUSS POINTS AND WEIGHTS
        x, w = gausslegendre(npoints)
        # RESCALING ON [binf,bsup]
        @. x = (bsup - binf)/2*x + (bsup + binf)/2
        @. w = (bsup - binf)/2*w
        # CACHE VARIABLES FOR FREE-ALLOCATION COMPUTATIONS
        wx = similar(x)
        shiftx = similar(wx)
        # EVALUATE POLYNOMIALS
        Qgenx = evaluate(basis.cache.prodMG, x)
        a = Int(sqrt(size(Qgenx,1)))
        Qgenxreshape = reshape(Qgenx, a, a, size(Qgenx,2))
        new{eltype(x)}(npoints,x,w,shiftx,wx, Qgenxreshape)
    end
end

