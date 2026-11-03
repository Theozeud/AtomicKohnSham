cd(@__DIR__)
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using AtomicKohnSham
using AtomicKohnSham: empty_channels, diagonalize_channels!, density_matrix!,
    assemble_hamiltonian!, QCFONCache
using CairoMakie
using LaTeXStrings
using LinearAlgebra
using Printf

# ODA (first order, energy line search) against QCFON (Newton on the density
# matrix with fractional occupations) on a few neutral atoms in LSDA, with the
# same discretization and tolerance. Both are measured with the same quantity,
# computed independently of the algorithm at every iteration: the
# self-consistency residual ‖D_out[D] - D‖ of the current density matrix D, where
# D_out[D] is built from the eigenvectors of H[D] with the current occupations.
# The energy error is measured against the QCFON energy. Produces:
#   - convergence_iterations.pdf: residual and energy error vs iteration;
#   - convergence_time.pdf: residual vs wall time (time spent in the monitoring
#     callback excluded, compilation excluded by a first run, fastest of 3 runs);
#   - comparison.txt: iterations, times, final residuals, energy differences.

# ===================== STEP 1 : PARAMETERS =====================#
atoms       = [(name = "C", Z = 6), (name = "Ne", Z = 10),
               (name = "Na", Z = 11), (name = "Cl", Z = 17)]

# DISCRETIZATION
Rmax        = 100
Nmesh       = 25
ordermax    = 10
lh          = 2
nh          = 4

# ALGORITHMS
scftol      = 1e-11
maxiter     = 500
oda   = ODA(tinit = 0.6, scftol = scftol, aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5))
qcfon = QCFON(scftol = scftol)

# ===================== STEP 2 : MONITORING =====================#
"""
    Monitor()

Callback recording, at every iteration, the self-consistency residual
‖D_out[D] - D‖ of the current iterate, its total energy, the wall time since the
first iteration (minus the time spent in this callback), and whether the
iteration is a Newton one.
"""
struct Monitor
    residual::Vector{Float64}
    energy::Vector{Float64}
    time::Vector{Float64}
    newton::Vector{Bool}
    clock::Vector{Float64}
end
Monitor() = Monitor(Float64[], Float64[], Float64[], Bool[], [0.0, 0.0])

"""
    true_residual(solver) -> Float64

‖D_out[D] - D‖ for the density matrix `D` and the occupations of the current
iterate, independently of the algorithm.
"""
function true_residual(solver::KSESolver)
    dis, model, c = solver.discretization, solver.model, solver.algcache
    nb, nspin = dis.Nₕ, dis.nspin
    D = reshape(copy(c.D), nb, nb, nspin)
    assemble_hamiltonian!(dis, model, c.D)
    channels = diagonalize_channels!(empty_channels(dis), dis, dis.ksham.H)
    for ch in channels
        nk = c.n[ch.l + 1, :, ch.σ]
        ch.occ = findall(>(0), nk)
        ch.nocc = nk[ch.occ]
    end
    Dout = zeros(nb, nb, nspin)
    density_matrix!(Dout, channels)
    norm(Dout .- D)
end

function AtomicKohnSham.callback!(m::Monitor, solver::KSESolver)
    tin = time()
    isempty(m.time) && (m.clock[1] = tin)
    push!(m.time, tin - m.clock[1] - m.clock[2])
    push!(m.residual, true_residual(solver))
    push!(m.energy, solver.energies.Etot)
    c = solver.algcache
    push!(m.newton, c isa QCFONCache && c.newton)
    m.clock[2] += time() - tin
    nothing
end

# ===================== STEP 3 : SOLVE =====================#
function problem(Z)
    ex = BuiltinFunctional(:lda_x; nspin = 2)
    ec = BuiltinFunctional(:lda_c_pw; nspin = 2)
    model = KSEModel(; Z = Z, N = Z, ex = ex, ec = ec)
    basis = P1IntLegendreBasis(expmesh(0, Rmax, Nmesh; s = 1.2); ordermax = ordermax)
    model, KSEDiscretization(basis, model; lh = lh, nh = nh)
end

"""
    run(Z, alg; nrep = 3)

Solve once to compile, then `nrep` times; return the solution and the monitor of
the fastest run (wall times on a shared machine fluctuate by a factor up to 2).
"""
function run(Z, alg; nrep = 3)
    model, dis = problem(Z)
    groundstate(model, dis, alg; maxiter = maxiter, callback = Monitor())
    best = nothing
    for _ in 1:nrep
        m = Monitor()
        sol = groundstate(model, dis, alg; maxiter = maxiter, callback = m)
        (best === nothing || m.time[end] < best[2].time[end]) && (best = (sol, m))
    end
    best
end

results = map(atoms) do a
    so, mo = run(a.Z, oda)
    sq, mq = run(a.Z, qcfon)
    (atom = a, oda = (sol = so, mon = mo), qcfon = (sol = sq, mon = mq))
end

# ===================== STEP 4 : COMPARISON.TXT =====================#
open("comparison.txt", "w") do io
    println(io, "Neutral atoms, LSDA (Slater + PW92), Rmax = $Rmax, $Nmesh cells, degree $ordermax, scftol = $scftol")
    println(io, "residual = ‖D_out[D] - D‖ at the last iterate; time excludes compilation and monitoring, fastest of 3 runs\n")
    @printf(io, "%-4s | %5s %8s %10s | %5s %6s %8s %10s | %10s\n", "atom", "ODA", "time(s)",
            "residual", "QCFON", "Newton", "time(s)", "residual", "E_ODA-E_QC")
    for r in results
        mo, mq = r.oda.mon, r.qcfon.mon
        @printf(io, "%-4s | %5d %8.2f %10.1e | %5d %6d %8.2f %10.1e | %10.1e\n", r.atom.name,
                length(mo.time), mo.time[end], mo.residual[end], length(mq.time),
                count(mq.newton), mq.time[end], mq.residual[end],
                r.oda.sol.energies.Etot - r.qcfon.sol.energies.Etot)
    end
end
print(read("comparison.txt", String))

# ===================== STEP 5 : PLOTS =====================#
floor_at(x) = max.(x, 1e-16)
colors = (oda = Makie.wong_colors()[1], qcfon = Makie.wong_colors()[2])

"""Residual and energy error of both algorithms for one atom, against `x(monitor)`."""
function panel!(ax, r, x, y)
    for (key, label) in ((:oda, "ODA"), (:qcfon, "QCFON"))
        m = getfield(r, key).mon
        X, Y = x(m), floor_at(y(m, r))
        lines!(ax, X, Y; color = colors[key], label = label)
        nw = findall(m.newton)
        key == :qcfon && scatter!(ax, X[nw], Y[nw]; color = colors[key], markersize = 7,
                                  label = "Newton iterations")
    end
end

residual(m, r) = m.residual
energy_error(m, r) = abs.(m.energy .- r.qcfon.sol.energies.Etot)

fig = Figure(size = (1000, 480))
for (j, r) in enumerate(results)
    ax1 = Axis(fig[1, j]; yscale = log10, title = r.atom.name,
               ylabel = j == 1 ? L"\Vert D_{\mathrm{out}}[D]-D\Vert" : "")
    panel!(ax1, r, m -> 1:length(m.time), residual)
    ax2 = Axis(fig[2, j]; yscale = log10, xlabel = "iteration",
               ylabel = j == 1 ? L"|E_k-E^*|" : "")
    panel!(ax2, r, m -> 1:length(m.time), energy_error)
end
Legend(fig[3, :], content(fig[1, 1]); orientation = :horizontal, merge = true, framevisible = false)
save("convergence_iterations.pdf", fig)

fig = Figure(size = (1000, 280))
for (j, r) in enumerate(results)
    ax = Axis(fig[1, j]; yscale = log10, title = r.atom.name, xlabel = "time (s)",
              ylabel = j == 1 ? L"\Vert D_{\mathrm{out}}[D]-D\Vert" : "")
    panel!(ax, r, m -> m.time, residual)
end
Legend(fig[2, :], content(fig[1, 1]); orientation = :horizontal, merge = true, framevisible = false)
save("convergence_time.pdf", fig)
