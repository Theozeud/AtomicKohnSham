# Cost of one application of the QCFON Newton operator (I - χ₀K, bordered) at a
# converged state, i.e. of one GMRES iteration, and of a whole QCFON run.
# Run: julia bench_operator.jl [label]
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", "examples"))
using AtomicKohnSham, LinearAlgebra, Printf
using AtomicKohnSham: newton_operator, prepare_newton_operator!

const LABEL = get(ARGS, 1, "run")

function solver_for(Z, N, nspin)
    ex = BuiltinFunctional(:lda_x; nspin = nspin)
    ec = BuiltinFunctional(:lda_c_pw; nspin = nspin)
    model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)
    basis = P1IntLegendreBasis(expmesh(0, 100, 20; s = 1.2); ordermax = 8)
    dis = KSEDiscretization(basis, model; lh = 2, nh = 4,
                            fem_integration_method = GaussLegendre(basis, 300))
    KSESolver(model, dis, QCFON(scftol = 1e-11); maxiter = 300)
end

for (name, Z, N, nspin) in (("Ne LDA", 10, 10, 1), ("Cl-0.5 LSDA", 17, 17.5, 2))
    solver = solver_for(Z, N, nspin)
    t = @elapsed AtomicKohnSham.solve!(solver)
    cache = solver.algcache
    nb = solver.discretization.Nₕ
    np = length(cache.partial) < 2 ? 0 : length(cache.partial) + 1
    x = randn(nb * nb * nspin + np)
    prepare_newton_operator!(cache, solver)
    newton_operator(x, cache, solver)
    top = @elapsed for _ in 1:10
        newton_operator(x, cache, solver)
    end
    @printf("%-12s [%s] run: %3d it %5.2fs | one operator application: %6.1f ms\n", name,
            LABEL, solver.niter, t, 100top)
end
