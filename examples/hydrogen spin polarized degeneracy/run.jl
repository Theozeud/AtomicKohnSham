cd(@__DIR__)
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using AtomicKohnSham
using CairoMakie
using Libxc

# Hydrogen, Z=N=1, spin-polarized (nspin=2), lh=0/nh=1 (1s only) -- the
# smallest possible case that isolates a single twofold degeneracy: the
# 1sUP/1sDOWN pair competing for 1 electron. Every model here has both spin
# channels seeing the *same* Hamiltonian at the very first iteration (a
# symmetric starting density), so `resolve_twofold_degeneracy!`
# (src/algorithms/aufbau/optimized.jl) is exercised on every single SCF step,
# not just the first.
#
# Two axes of comparison:
#   - functional: `hartree_only` (ex=ec=NoFunctional) vs `lda` (Slater
#     exchange + Perdew-Wang correlation). For `hartree_only`, nothing in the
#     energy functional distinguishes spin at all (Hartree only sees the
#     total density) -- ε(1sUP) = ε(1sDOWN) EXACTLY, every iteration, not
#     just at start. This is the purest possible reproduction of the "flat
#     energy w.r.t. occupation split" pathology. For `lda`, Slater exchange
#     genuinely prefers full spin polarization (a single electron split
#     50/50 across both spin channels has *higher* exchange energy than the
#     same electron fully in one channel -- a real, not roundoff-sized,
#     energy difference), so there IS a correct answer to find.
#   - guard on/off: `flat_tol` on `OptimizedAufbau` controls whether
#     `resolve_twofold_degeneracy!` falls back to a symmetric occupation
#     split when the sampled energy is flat (see docstring), or always
#     trusts `line_search_energy`'s minimizer -- `flat_tol=0` reproduces the
#     old, unguarded behavior for comparison.
#   - `hartree_only_smeared` swaps `OptimizedAufbau` for `SmearedAufbau`
#     entirely, sidestepping the discrete degeneracy-resolution machinery.
#
# See description.md for what was actually observed running this.

# ===================== STEP 1 : PARAMETERS =====================#
Z, N = 1, 1
Rmax        = 1000
Nmesh       = 30
ordermax    = 10
lh          = 0
nh          = 1
scftol      = 1e-10
maxiter     = 300   # generous: the guard-off cases are exactly what might not converge cleanly

# --- per-iteration occupation tracker, so we can see n(1sUP)/n(1sDOWN) settle
#     (or fail to) rather than only the final answer ---
struct OccTracker
    up::Vector{Float64}
    down::Vector{Float64}
end
OccTracker() = OccTracker(Float64[], Float64[])
AtomicKohnSham.callback!(cb::OccTracker, solver::KSESolver) = begin
    push!(cb.up, Float64(solver.algcache.n[1, 1, 1]))
    push!(cb.down, Float64(solver.algcache.n[1, 1, 2]))
end

# ===================== STEP 2 : BUILD & SOLVE =====================#
cases = [
    (key = "hartree_only_guard_on",  name = "Hartree only, guard ON",
     ex = NoFunctional(2), ec = NoFunctional(2),
     aufbau = OptimizedAufbau(max_degen = 2)),
    (key = "hartree_only_guard_off", name = "Hartree only, guard OFF",
     ex = NoFunctional(2), ec = NoFunctional(2),
     aufbau = OptimizedAufbau(max_degen = 2, flat_tol = 0)),
    (key = "lda_guard_on",           name = "LDA, guard ON",
     ex = Functional(:lda_x, n_spin = 2), ec = Functional(:lda_c_pw, n_spin = 2),
     aufbau = OptimizedAufbau(max_degen = 2)),
    (key = "lda_guard_off",          name = "LDA, guard OFF",
     ex = Functional(:lda_x, n_spin = 2), ec = Functional(:lda_c_pw, n_spin = 2),
     aufbau = OptimizedAufbau(max_degen = 2, flat_tol = 0)),
    (key = "hartree_only_smeared",   name = "Hartree only, SmearedAufbau",
     ex = NoFunctional(2), ec = NoFunctional(2),
     aufbau = SmearedAufbau(Temp = 0.01)),
]

sols = KSESolution[]
trackers = OccTracker[]
for c in cases
    model = KSEModel(; Z = Z, N = N, hartree = 1, ex = c.ex, ec = c.ec)

    mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
    dis = KSEDiscretization(basis, model; lh = lh, nh = nh)

    alg = ODA(tinit = 0.6, aufbau = c.aufbau, scftol = scftol)
    tracker = OccTracker()

    sol = open(joinpath(c.key, "log.txt"), "w") do io
        write_log_header(io, model, alg, dis)
        groundstate(model, dis, alg; maxiter = maxiter,
            callback = CallbackSet((LogFileCallback(io), tracker)))
    end

    println(c.name, ": ", sol.success, " in ", sol.niter, " iterations",
        " (final stopping criterion = ", sol.stopping_criteria, ")")
    write_report(sol, joinpath(c.key, "results.txt"))
    push!(sols, sol)
    push!(trackers, tracker)
end

# ===================== STEP 3 : COMPARISON.TXT =====================#
const COLW = maximum(length(c.name) for c in cases) + 2
_col(s, w = COLW) = rpad(s, w)
_round0(v; digits = 6) = (r = round(v, digits = digits); iszero(r) ? zero(r) : r)

open("comparison.txt", "w") do io
    println(io, "Hydrogen (Z=1, N=1), spin-polarized -- 1sUP/1sDOWN degeneracy resolution")
    println(io, "="^78)
    println(io)

    println(io, "SCF OUTCOME")
    println(io, "-"^78)
    for (c, sol) in zip(cases, sols)
        println(io, _col(c.name), sol.success, " in ", sol.niter,
            " iterations (stopping criterion = ", sol.stopping_criteria, ")")
    end
    println(io)

    println(io, "ENERGIES (Ha)")
    println(io, "-"^78)
    println(io, _col("", 12), join(_col(c.name) for c in cases))
    for f in (:Ekin, :Ecou, :Ehar, :Eexc, :Etot)
        vals = [getfield(sol.energies, f) for sol in sols]
        println(io, _col(string(f), 12), join(_col(string(_round0(v))) for v in vals))
    end
    println(io)

    println(io, "FINAL OCCUPATION (n(1sUP), n(1sDOWN))")
    println(io, "-"^78)
    for (c, tr) in zip(cases, trackers)
        println(io, _col(c.name), "n_up = ", round(tr.up[end]; digits = 6),
            "   n_down = ", round(tr.down[end]; digits = 6))
    end
    println(io)

    println(io, "OCCUPATION STABILITY (spread of n(1sUP) over the LAST 20 iterations)")
    println(io, "-"^78)
    for (c, tr) in zip(cases, trackers)
        tail = tr.up[max(1, end - 19):end]
        println(io, _col(c.name), "max-min(n_up, last 20 iters) = ",
            round(maximum(tail) - minimum(tail); digits = 6))
    end
end

# ===================== STEP 4 : COMPARISON PLOTS =====================#
const FS_LABEL, FS_TICK, FS_LEGEND, LW = 26, 18, 16, 4
const PALETTE = [:orange, :seagreen, :dodgerblue, :firebrick, :royalblue]

_axis(fig_pos; kwargs...) = Axis(fig_pos; xlabelsize = FS_LABEL, ylabelsize = FS_LABEL,
    xticklabelsize = FS_TICK, yticklabelsize = FS_TICK, kwargs...)

# --- stopping criterion vs iteration: clean decay vs stalled/oscillating ---
fig = Figure(size = (900, 650))
ax = _axis(fig[1, 1]; xlabel = "Iteration", ylabel = "Stopping criterion", yscale = log10)
for (i, (c, sol)) in enumerate(zip(cases, sols))
    sc = max.(sol.logbook.stopping_criteria, eps(Float64))
    lines!(ax, eachindex(sc), sc; linewidth = LW, label = c.name, color = PALETTE[i])
end
hlines!(ax, [scftol]; linewidth = 2, linestyle = :dash, color = :black, label = "scftol")
axislegend(ax; labelsize = FS_LEGEND, position = :rt)
save("convergence_comparison.pdf", fig)

# --- n(1sUP) vs iteration: settles to a fixed point, or keeps moving ---
fig = Figure(size = (900, 650))
ax = _axis(fig[1, 1]; xlabel = "Iteration", ylabel = "n(1sUP)")
for (i, (c, tr)) in enumerate(zip(cases, trackers))
    lines!(ax, eachindex(tr.up), tr.up; linewidth = LW, label = c.name, color = PALETTE[i])
end
axislegend(ax; labelsize = FS_LEGEND, position = :rb)
save("occupation_comparison.pdf", fig)
