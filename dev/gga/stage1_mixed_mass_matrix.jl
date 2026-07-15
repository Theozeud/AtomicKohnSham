# Stage 1 (see PLAN.md): validates the new Mmix[w] = ∫ w(r) φᵢ'(r)φⱼ(r) dr
# FEM primitive against direct numerical quadrature (QuadGK), independent of
# the FEM machinery being tested. Cross-checks a handful of (i,j) entries for
# an arbitrary smooth weight w(r), on a non-uniform mesh (to make sure the
# per-cell Jacobian/chain-rule bookkeeping -- which this session's derivation
# claims cancels exactly for this one-derivative building block -- is right).
#
# Run: julia --project=. dev/gga/stage1_mixed_mass_matrix.jl

using AtomicKohnSham
using QuadGK
using Test

Rmax = 12.0
mesh = expmesh(0, Rmax, 6; s = 1.3)
basis = P1IntLegendreBasis(mesh; ordermax = 4)

w(r) = exp(-0.3r) * (1 + 0.1r^2)   # arbitrary smooth test weight

weight = AtomicKohnSham.FunWeight(w; is_vectorized = false, is_inplace = false)
Mmix = mixed_mass_matrix(basis; weight = weight)

@testset "Mmix[w] vs direct numerical quadrature" begin
    n = length(basis)
    # sample a handful of (i,j) pairs, not the full n^2 (QuadGK per-pair is slow)
    pairs = [(i, j) for i in 1:3:n for j in 1:3:n]
    for (i, j) in pairs
        f(r) = begin
            h = 1e-6
            dphi_i = (basis(i, r + h) - basis(i, r - h)) / (2h)
            dphi_i * basis(j, r) * w(r)
        end
        val, _ = quadgk(f, 0.0, Rmax; rtol = 1e-10)
        @test isapprox(Mmix[i, j], val; atol = 1e-7, rtol = 1e-6)
    end
end
