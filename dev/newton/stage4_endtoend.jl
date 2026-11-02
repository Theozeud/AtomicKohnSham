# Stage 4: QCFON against ODA (same discretization, same tolerance), with the
# residual history, on atoms and ions with and without a fractional orbital.
# Run: julia stage4_endtoend.jl [switch_tol]
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", "examples"))
using AtomicKohnSham, Printf

"""
Records the stopping criterion and the phase of each iteration: an iteration is a
Newton one if Newton had already started before it (the ODA iteration that
triggers the switch is an ODA one).
"""
struct History
    crit::Vector{Float64}
    newton::Vector{Bool}
    started::Base.RefValue{Bool}
end
History() = History(Float64[], Bool[], Ref(false))
function AtomicKohnSham.callback!(h::History, solver::KSESolver)
    push!(h.crit, Float64(solver.stopping_criteria))
    push!(h.newton, h.started[])
    h.started[] = solver.algcache isa AtomicKohnSham.QCFONCache && solver.algcache.newton
end

function run(Z, N, nspin, alg; lh = 2, nh = 4)
    ex = BuiltinFunctional(:lda_x; nspin = nspin)
    ec = BuiltinFunctional(:lda_c_pw; nspin = nspin)
    model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)
    basis = P1IntLegendreBasis(expmesh(0, 60, 15; s = 1.2); ordermax = 8)
    dis = KSEDiscretization(basis, model; lh = lh, nh = nh,
                            fem_integration_method = GaussLegendre(basis, 300))
    h = History()
    t = @elapsed sol = groundstate(model, dis, alg; maxiter = 500, callback = h)
    sol, h, t
end

scftol = 1e-11
oda = ODA(tinit = 0.6, scftol = scftol, aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5))
qc = QCFON(scftol = scftol, switch_tol = parse(Float64, get(ARGS, 1, "1e-2")))
for (name, Z, N, nspin) in (("H LSDA", 1, 1, 2), ("Ne LDA", 10, 10, 1),
                            ("C LSDA", 6, 6, 2), ("Cl- 0.5 LSDA", 17, 17.5, 2))
    so, ho, to = run(Z, N, nspin, oda)
    sq, hq, tq = run(Z, N, nspin, qc)
    @printf("%-13s ODA: %3d it %5.1fs E=%.12f | QCFON: %3d it (%d Newton) %5.1fs E=%.12f | ΔE=%.1e\n",
            name, so.niter, to, so.energies.Etot, sq.niter, count(hq.newton), tq,
            sq.energies.Etot, sq.energies.Etot - so.energies.Etot)
    println("    Newton residuals: ", join([@sprintf("%.1e", c) for c in hq.crit[hq.newton]], "  "))
end
