cd(@__DIR__)
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using AtomicKohnSham
using CairoMakie
using Libxc

# ===================== STEP 1 : PARAMETERS =====================#
# MODEL
Z = 11
N = 11
ex = Functional(:lda_x; n_spin = 1)
ec = NoFunctional(1)

# DISCRETIZATION : P1IntLegendreBasis + expmesh
Rmax        = 500
Nmesh       = 20
ordermax    = 10
lh          = 2

# ALGORITHM : ODA + Optimized Aufbau
maxiter     = 100
degen_tol   = 1e-1
scftol      = 1e-9

# PLOT
Rmax_plot   = 300
N_plot      = 2000

# ===================== STEP 2 : BUILD & SOLVE =====================#
model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)

mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
dis = KSEDiscretization(basis, model; lh = lh, nh = 10,
    fem_integration_method = GaussLegendre(basis, 2000))

aufbau = OptimizedAufbau(max_degen = 2, tol = degen_tol)
alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = scftol)

# log.txt is written incrementally during the SCF loop through a
# LogFileCallback, so it stays useful for investigation even if the
# computation is interrupted or throws.
sol = open("log.txt", "w") do io
    write_log_header(io, model, alg, dis)
    groundstate(model, dis, alg; maxiter = maxiter,
        callback = CallbackSet((LogFileCallback(io),)))
end

println(sol)

# ===================== STEP 3 : RESULTS.TXT =====================#
write_report(sol, "results.txt")

# ===================== STEP 4 : PLOTS =====================#
X = AtomicKohnSham.exprange(1e-3, Rmax_plot, N_plot; s = 1.5)

save("density.pdf", plot_density(sol, X))
save("orbitals.pdf", plot_orbitals(sol, X))
save("potentials.pdf", plot_potentials(sol, X; l = 0))
save("convergence.pdf", plot_convergence(sol))
save("energy_breakdown.pdf", plot_energy_breakdown(sol))
