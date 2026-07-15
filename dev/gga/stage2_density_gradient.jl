# Stage 2 (see PLAN.md): validates eval_density_gradient!/eval_density_gradient
# (the general per-point evaluator) against finite differences of the
# existing eval_density!, and cross-checks the "optimized" per-cell-quadrature
# version (optimized_eval_density_gradient!, used inside assemble_exc!) agrees
# with the general one -- both on a REAL converged density matrix from a
# Hydrogen SCF run (Slater exchange, so a genuinely non-trivial density).
#
# Run: julia --project=. dev/gga/stage2_density_gradient.jl

using AtomicKohnSham
using Test

Z, N = 1.0, 1.0
Rmax, Nmesh, ordermax = 30.0, 20, 8
ex = BuiltinFunctional(:lda_x; nspin = 1)
ec = NoFunctional(1)
model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)
mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
dis = KSEDiscretization(basis, model; lh = 0, nh = 1)
aufbau = OptimizedAufbau(max_degen = 1)
alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-10)
sol = groundstate(model, dis, alg; maxiter = 100)
D = sol.D

X = collect(range(0.05, Rmax - 0.5; length = 25))

@testset "eval_density_gradient! vs finite difference of eval_density" begin
    ρ, dρ = AtomicKohnSham.eval_density_gradient(dis, D, X)
    h = 1e-5
    for (k, x) in enumerate(X)
        ρm = AtomicKohnSham.eval_density(dis, D, x - h)
        ρp = AtomicKohnSham.eval_density(dis, D, x + h)
        fd = (ρp - ρm) / (2h)
        @test isapprox(ρ[k], AtomicKohnSham.eval_density(dis, D, x); rtol = 1e-12)
        @test isapprox(dρ[k], fd; atol = 1e-6, rtol = 1e-5)
    end
end

@testset "optimized_eval_density_gradient! agrees with the general one (per cell)" begin
    # optimized_eval_density!/optimized_eval_density_gradient! assume X is the
    # *current cell's* shifted quadrature nodes (Qgenx/Qmixedgenx are indexed
    # by quadra.x's position, not by X's values) -- exactly how assemble_exc!
    # calls them via fill_local_matrix!'s shiftx. Reconstruct that for one cell.
    quadra = dis.fem_integration_method
    cellidx = 3
    invϕ = basis.invshifts[cellidx]
    shiftx = invϕ[1] .* quadra.x .+ invϕ[2]
    ρ_gen, dρ_gen = AtomicKohnSham.eval_density_gradient(dis, D, shiftx)
    ρ_opt = zeros(length(shiftx))
    AtomicKohnSham.optimized_eval_density!(ρ_opt, dis, D, shiftx)
    dρ_opt = zeros(length(shiftx))
    AtomicKohnSham.optimized_eval_density_gradient!(dρ_opt, dis, D, ρ_opt, shiftx)
    @test isapprox(ρ_opt, ρ_gen; rtol = 1e-10)
    @test isapprox(dρ_opt, dρ_gen; rtol = 1e-8)
end

println("Sample: x=$(X[5]) ρ=$(AtomicKohnSham.eval_density(dis,D,X[5]))")
