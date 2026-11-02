# Stage 5: iron in reduced Hartree-Fock (Z = N = 26, no exchange-correlation,
# spin-unpolarized). The Fermi level is a degenerate 5s/3d level very close to 0
# (Cancès-Mourad, CAMCoS 2018: n(3d) = 5 x 1.1957, eps_F ~ -1e-5 Ha), so it tests
# QCFON with two fractional orbitals in two channels. nh must include the 5s
# (nh >= 5); with nh = 4 the 5s cannot be occupied and the 3d is pushed above 0.
# Run: julia stage5_iron_rhf.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", "examples"))
using AtomicKohnSham, Printf

"""Records the stopping criterion and whether Newton had started before the iteration."""
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

function run(alg; maxiter = 400)
    model = RHF(Z = 26, N = 26)
    basis = P1IntLegendreBasis(expmesh(0, 240, 30; s = 1.2); ordermax = 8)
    dis = KSEDiscretization(basis, model; lh = 2, nh = 8,
                            fem_integration_method = GaussLegendre(basis, 300))
    h = History()
    t = @elapsed sol = groundstate(model, dis, alg; maxiter = maxiter, callback = h)
    sol, h, t
end

scftol = 1e-11
results = [("ODA", run(ODA(tinit = 0.6, scftol = scftol,
                           aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5)))),
           ("QCFON", run(QCFON(scftol = scftol)))]
for (name, (sol, h, t)) in results
    @printf("%-6s %s  %3d it (%d Newton)  %6.1fs  E = %.12f  final crit = %.1e\n", name,
            sol.success, sol.niter, count(h.newton), t, sol.energies.Etot, h.crit[end])
    for (lab, ϵ, n) in sol.occupied[max(1, end - 2):end]
        @printf("        %-4s ϵ = % .10f  n = %.10f\n", lab, ϵ, n)
    end
    any(h.newton) && println("        Newton residuals: ",
        join([@sprintf("%.1e", c) for c in h.crit[h.newton]], "  "))
    println("        last 6 criteria: ", join([@sprintf("%.1e", c) for c in h.crit[max(1, end-5):end]], "  "))
end
