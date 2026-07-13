cd(@__DIR__)
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using AtomicKohnSham
using CairoMakie
using DoubleFloats

# This is examples/basic_example redone with all FEM/SCF arithmetic in
# Double64 (~32 significant digits) instead of Float64, with the same
# physical/discretization/algorithm parameters EXCEPT the exchange
# functional -- see description.md for why: exchange-correlation energy is
# evaluated by Gauss-Legendre quadrature (GaussLegendre/FastGaussQuadrature),
# which is hardcoded to Float64 with no generic-precision variant, so it
# silently caps Etot at Float64 precision regardless of everything else
# being Double64. Hartree=1 is kept (the Hartree term is exact polynomial
# integration, not quadrature, so it doesn't have this problem) but there is
# no exchange-correlation functional here, unlike basic_example's Slater
# exchange -- this trades "same physics" for "actually demonstrates Double64
# precision", which running with Slater exchange included cannot do until
# GaussLegendre gets a generic-precision node/weight generator.
#
# Z, N are also kept as plain Int rather than Double64: an unrelated
# internal Int-only signature on the twofold-degeneracy resolution path
# would otherwise reject them (see description.md).

# ===================== STEP 1 : PARAMETERS =====================#
# MODEL
Z = 11
N = 11
ex = NoFunctional(1)
ec = NoFunctional(1)

# DISCRETIZATION : P1IntLegendreBasis + expmesh
Rmax        = 500
Nmesh       = 20
ordermax    = 10
lh          = 2

# ALGORITHM : ODA + Optimized Aufbau
maxiter     = 100
degen_tol   = Double64(1e-2)
scftol      = Double64(1e-20)

# PLOT
Rmax_plot   = 300
N_plot      = 2000

# ===================== STEP 2 : BUILD & SOLVE =====================#
model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)

mesh = expmesh(0, Rmax, Nmesh; T = Double64, s = 1.2)
basis = P1IntLegendreBasis(mesh, Double64; ordermax = ordermax)
# GaussLegendre is still built here (matching basic_example structurally),
# but with no exchange-correlation functional it's never actually touched by
# the computation -- Ekin/Ecou/Ehar all go through exact integration.
dis = KSEDiscretization(basis, model; lh = lh, nh = 10,
    fem_integration_method = GaussLegendre(basis, 2000))

aufbau = OptimizedAufbau(max_degen = 2, tol = degen_tol)
alg = ODA(tinit = Double64(0.6), aufbau = aufbau, scftol = scftol)

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
# Plotted in Float64: CairoMakie renders to pixels, so the extra precision
# of X/the solution buys nothing visually. eval_density/eval_orbital/etc.
# promote X to the solution's Double64 internally regardless.
X = AtomicKohnSham.exprange(1e-3, Rmax_plot, N_plot; s = 1.5)

save("density.pdf", plot_density(sol, X))
save("orbitals.pdf", plot_orbitals(sol, X))
save("potentials.pdf", plot_potentials(sol, X; l = 0))
save("convergence.pdf", plot_convergence(sol))
save("energy_breakdown.pdf", plot_energy_breakdown(sol))
