# Benchmark of the quadrature of the exchange-correlation energy (LDA/LSDA).
#
# 1. Integration error alone: Exc of one fixed, converged density matrix,
#    computed with GaussLegendre(basis, n) for several n, against n = 4000.
# 2. Full SCF (ODA): iterations, time per iteration, total energy and Fermi level
#    against the n = 2000 run.
#
# Usage: julia benchmark/exc_energy_quadrature.jl [label]
# (the label, e.g. "before"/"after", only tags the output file).
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "examples"))
using AtomicKohnSham, Printf
using AtomicKohnSham: compute_exc_energy

const LABEL = get(ARGS, 1, "run")
const NQUADS = (2000, 500, 200, 100, 50, 30, 20)
const CASES = (("Ne LDA", 10, 10, 1), ("Cl LSDA", 17, 17, 2), ("Ar-0.45 LSDA", 18, 18.45, 2))

function problem(Z, N, nspin, nquad)
    ex = BuiltinFunctional(:lda_x; nspin = nspin)
    ec = BuiltinFunctional(:lda_c_pw; nspin = nspin)
    model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)
    basis = P1IntLegendreBasis(expmesh(0, 100, 20; s = 1.2); ordermax = 8)
    dis = KSEDiscretization(basis, model; lh = 2, nh = 4,
                            fem_integration_method = GaussLegendre(basis, nquad))
    model, dis
end

const ALG = ODA(tinit = 0.6, scftol = 1e-11, aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5))

out = IOBuffer()
for (name, Z, N, nspin) in CASES
    model, dis = problem(Z, N, nspin, 2000)
    ref = groundstate(model, dis, ALG; maxiter = 500)
    _, disref = problem(Z, N, nspin, 4000)
    Eref = compute_exc_energy(disref, model, ref.D)
    println(out, "== $name: integration error of Exc on a fixed density matrix (vs n = 4000)")
    for nq in NQUADS
        _, d = problem(Z, N, nspin, nq)
        compute_exc_energy(d, model, ref.D)
        t = @elapsed for _ in 1:5; compute_exc_energy(d, model, ref.D); end
        @printf(out, "   n = %4d   |ΔExc| = %.1e   time = %6.2f ms\n", nq,
                abs(compute_exc_energy(d, model, ref.D) - Eref), 1000t / 5)
    end
    println(out, "== $name: full SCF (vs n = 2000)")
    ϵFref = last(ref.occupied)[2]
    for nq in NQUADS
        m, d = problem(Z, N, nspin, nq)
        t = @elapsed sol = groundstate(m, d, ALG; maxiter = 500)
        @printf(out, "   n = %4d   %-8s %3d it  %6.1f ms/it  ΔE = % .1e   ΔϵF = % .1e\n", nq,
                sol.success, sol.niter, 1000t / sol.niter, sol.energies.Etot - ref.energies.Etot,
                last(sol.occupied)[2] - ϵFref)
    end
end
s = String(take!(out))
print(s)
mkpath(joinpath(@__DIR__, "results"))
write(joinpath(@__DIR__, "results", "exc_energy_quadrature_$LABEL.txt"), s)
