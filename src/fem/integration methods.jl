abstract type FEMIntegrationMethod end

struct NoSelectedMethod <: FEMIntegrationMethod end

struct ExactIntegration <: FEMIntegrationMethod end

abstract type QuadratureIntegration <: FEMIntegrationMethod end

"""
    GaussLegendre(basis::FEMBasis, npoints::Int = 1000)

Gauss–Legendre quadrature method for FEM matrix assembly and energy
integrals, with `npoints` points per mesh cell.

Precomputes and caches the quadrature nodes/weights (rescaled to each
reference cell and to the full radial domain `[0, Rmax]`), together with the
FEM generator polynomials evaluated at the nodes, so that repeated
assemblies during the SCF loop avoid recomputing them.

The nodes/weights are generated at `eltype(basis)` precision: for `Float64`
(the common case) via `FastGaussQuadrature.gausslegendre`, an O(`npoints`)
algorithm that is Float64-only; for any other type (e.g. `Double64`,
`BigFloat`) via `QuadGK.gauss`, a generic-precision O(`npoints`²)
Golub–Welsch implementation. The latter is markedly slower at equal
`npoints` (benchmarked ~1s for `Double64` at `npoints=2000`), but this
constructor only runs once per discretization, not per SCF iteration, so
that one-time cost is worth paying to get genuinely `T`-accurate nodes
rather than nodes that only ever had Float64 accuracy to begin with.
"""
struct GaussLegendre{T <: Real} <: FEMIntegrationMethod
    npoints::Int
    # FOR FEM MATRIX COMPUTATIONS
    x::Vector{T}
    w::Vector{T}
    Pgenx::Matrix{T}
    Qgenx::Matrix{T}
    # FOR ENERGY COMPUTATIONS
    y::Vector{T}
    wy::Vector{T}
    # CACHE VARIABLES
    shiftx::Vector{T}
    fx::Vector{T}
    fy::Vector{T}
    fx2::Matrix{T}
    function GaussLegendre(basis::FEMBasis, npoints::Int = 1000)
        @unpack binf, bsup = basis.generators
        T = eltype(basis)
        # GENERATION OF GAUSS POINTS AND WEIGHTS ON [-1,1]
        x, w = T === Float64 ? gausslegendre(npoints) : quadgk_gauss(T, npoints)
        # RESCALING ON [0,Rmax]
        Rmax = last(basis.mesh)
        y = Rmax/2 .* x .+ Rmax/2
        wy = Rmax/2 .* w
        # RESCALING ON [binf,bsup]
        @. x = (bsup - binf)/2*x + (bsup + binf)/2
        @. w = (bsup - binf)/2*w
        # CACHE VARIABLES FOR FREE-ALLOCATION COMPUTATIONS
        fx = similar(x)
        fy = similar(fx)
        shiftx = similar(fx)
        fx2 = similar(fx, 2, length(fx))
        # EVALUATE POLYNOMIALS
        Pgenx = evaluate(basis.generators.polynomials,x)
        Qgenx = evaluate(basis.cache.prodMG, x)
        #a = Int(sqrt(size(Qgenx, 1)))
        #Qgenxreshape = reshape(Qgenx, a, a, size(Qgenx, 2))
        new{T}(npoints, x, w, Pgenx, Qgenx, y, wy, shiftx, fx, fy, fx2)
    end
end
