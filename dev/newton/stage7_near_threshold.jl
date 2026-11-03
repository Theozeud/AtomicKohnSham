# Stage 7: QCFON near the ionization threshold of chlorine (LSDA, Slater + PW92,
# Nc = 17.9864313 on this discretization: Rmax = 200, 30 cells, degree 10, the
# setup of dev/ionization_root). As N → Nc the Fermi level goes to 0 and box
# states come close to it, so the sum-over-states coefficients n_i/(ϵ_i - ϵ_a) grow
# and the Newton system may become ill-conditioned. For each N: Newton residuals,
# GMRES applications of the operator per Newton step, and timings against ODA.
# Run: julia stage7_near_threshold.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", "examples"))
using AtomicKohnSham, Printf, Logging

"""Collects the number of operator applications of each GMRES solve from the `@debug` messages of QCFON."""
struct GMRESLog <: AbstractLogger
    numops::Vector{Int}
end
Logging.min_enabled_level(::GMRESLog) = Logging.Debug
Logging.shouldlog(::GMRESLog, level, _module, args...) = _module === AtomicKohnSham
Logging.catch_exceptions(::GMRESLog) = false
function Logging.handle_message(l::GMRESLog, level, message, args...; kwargs...)
    message == "QCFON Newton step" && push!(l.numops, kwargs[Symbol("info.numops")])
    nothing
end

"""Newton residuals (stopping criterion of the Newton iterations)."""
struct Residuals
    crit::Vector{Float64}
end
function AtomicKohnSham.callback!(r::Residuals, solver::KSESolver)
    c = solver.algcache
    c isa AtomicKohnSham.QCFONCache && c.newton && push!(r.crit, solver.stopping_criteria)
end

function problem(N)
    ex = BuiltinFunctional(:lda_x; nspin = 2)
    ec = BuiltinFunctional(:lda_c_pw; nspin = 2)
    model = KSEModel(; Z = 17, N = N, ex = ex, ec = ec)
    basis = P1IntLegendreBasis(expmesh(0, 200, 30; s = 1.2); ordermax = 10)
    model, KSEDiscretization(basis, model; lh = 2, nh = 4)
end

scftol = 1e-11
oda = ODA(tinit = 0.6, scftol = scftol, aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5))
qc = QCFON(scftol = scftol)
let (m, d) = problem(17.5)
    groundstate(m, d, oda; maxiter = 1000)
    groundstate(m, d, qc; maxiter = 1000, callback = Residuals(Float64[]))
end
for N in (17.5, 17.9, 17.98, 17.9864, 17.98643, 17.9865, 17.99)
    model, dis = problem(N)
    to = @elapsed so = groundstate(model, dis, oda; maxiter = 1000)
    r, lg = Residuals(Float64[]), GMRESLog(Int[])
    tq = @elapsed sq = with_logger(lg) do
        groundstate(model, dis, qc; maxiter = 1000, callback = r)
    end
    ϵF = maximum(sq.ϵ[sq.n .> 0])
    @printf("N=%-9.5f ϵF=%+.3e | ODA %4d it %5.2fs | QCFON %3d it (%d Newton) %5.2fs | ΔE=%.1e\n",
            N, ϵF, so.niter, to, sq.niter, length(r.crit), tq,
            sq.energies.Etot - so.energies.Etot)
    println("    Newton residuals: ", join([@sprintf("%.1e", c) for c in r.crit], " "))
    println("    GMRES operator applications: ", join(lg.numops, " "))
end
