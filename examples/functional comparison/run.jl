cd(@__DIR__)
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using AtomicKohnSham
using CairoMakie
using Libxc

# Same atom (spin-unpolarized oxygen, Z=N=8), same discretization/algorithm,
# five different electron-electron interaction models -- to see how much of
# the total energy, the orbital spectrum, and the density comes from each
# physical term:
#   1. no_functional: Hartree term ON, no XC (mean-field electron-electron
#                     repulsion only)
#   2. lda_x_only   : Hartree ON + Slater exchange only, no correlation
#                     (BuiltinFunctional)
#   3. lda          : Hartree ON + Slater exchange + Perdew-Wang correlation
#                     (BuiltinFunctional, pure Julia, LDA)
#   4. lda_libxc    : same functionals as `lda` (Slater exchange + PW92
#                     correlation), but evaluated through Libxc instead of
#                     the BuiltinFunctional implementation -- a cross-check
#                     of the pure-Julia LDA against the reference C library
#   5. pbe          : Hartree ON + PBE exchange-correlation (Libxc, GGA)
# Each case gets its own log.txt/results.txt subfolder (same convention as
# examples/basic_example_doublefloat), and the top-level folder collects the
# side-by-side comparison: comparison.txt plus three overlay plots.
#
# groundstate() always names a solution from (Z, N) alone (here: "Oxygen"
# for every case), so that name can't be used to tell the cases apart -- the
# `cases` NamedTuples below carry the display names instead, used directly
# everywhere a label is needed (comparison.txt, the plots).

# ===================== STEP 1 : PARAMETERS =====================#
# MODEL
Z = 8
N = 8

# DISCRETIZATION : P1IntLegendreBasis + expmesh (same for every case, so the
# comparison isn't confounded by discretization error). lh=2/nh=5 is the
# same headroom used for scandium in test/scandium.jl -- generous for an
# atom whose ground state only fills up to 2p.
Rmax        = 1000
Nmesh       = 30
ordermax    = 10
lh          = 2
nh          = 5

# ALGORITHM : ODA + Optimized Aufbau (same for every case)
maxiter     = 100
scftol      = 1e-10

# PLOT
X = AtomicKohnSham.exprange(1e-3, 15, 2000; s = 1.5)

# ===================== STEP 2 : BUILD & SOLVE =====================#
cases = [
    (key = "no_functional", name = "Hartree only (no XC)",
     hartree = 1, ex = NoFunctional(1), ec = NoFunctional(1)),
    (key = "lda_x_only",    name = "LDA exchange only (Slater)",
     hartree = 1, ex = BuiltinFunctional(:lda_x; nspin = 1), ec = NoFunctional(1)),
    (key = "lda",           name = "LDA (Slater + Perdew-Wang)",
     hartree = 1, ex = BuiltinFunctional(:lda_x; nspin = 1), ec = BuiltinFunctional(:lda_c_pw; nspin = 1)),
    (key = "lda_libxc",     name = "LDA via Libxc (Slater + Perdew-Wang)",
     hartree = 1, ex = Functional(:lda_x, n_spin = 1), ec = Functional(:lda_c_pw, n_spin = 1)),
    (key = "pbe",           name = "PBE (GGA)",
     hartree = 1, ex = Functional(:gga_x_pbe, n_spin = 1), ec = Functional(:gga_c_pbe, n_spin = 1)),
]

sols = KSESolution[]
for c in cases
    model = KSEModel(; Z = Z, N = N, hartree = c.hartree, ex = c.ex, ec = c.ec)

    mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
    dis = KSEDiscretization(basis, model; lh = lh, nh = nh)

    aufbau = OptimizedAufbau(max_degen = 2)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = scftol)

    sol = open(joinpath(c.key, "log.txt"), "w") do io
        write_log_header(io, model, alg, dis)
        groundstate(model, dis, alg; maxiter = maxiter,
            callback = CallbackSet((LogFileCallback(io),)))
    end

    println(c.name, ": ", sol.success, " in ", sol.niter, " iterations")
    write_report(sol, joinpath(c.key, "results.txt"))
    push!(sols, sol)
end

# ===================== STEP 3 : COMPARISON.TXT =====================#
_round0(v; digits = 6) = (r = round(v, digits = digits); iszero(r) ? zero(r) : r)

const COLW = maximum(length(c.name) for c in cases) + 2
_col(s, w = COLW) = rpad(s, w)

open("comparison.txt", "w") do io
    println(io, "Oxygen (Z=$Z, N=$N), spin-unpolarized -- functional comparison")
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

    println(io, "OCCUPIED ORBITALS (orbital energy in Ha / occupation)")
    println(io, "-"^78)
    shells = sort(unique(occ[1] for sol in sols for occ in sol.occupied),
        by = s -> AtomicKohnSham.parse_shell(s)[1:2])
    println(io, _col("shell", 8), join(_col(c.name) for c in cases))
    for shell in shells
        row = _col(shell, 8)
        for sol in sols
            m = findfirst(occ -> occ[1] == shell, sol.occupied)
            cell = m === nothing ? "--" :
                "$(_round0(sol.occupied[m][2]; digits = 4)) / $(sol.occupied[m][3])"
            row *= _col(cell)
        end
        println(io, row)
    end
end

# ===================== STEP 4 : COMPARISON PLOTS =====================#
const FS_LABEL, FS_TICK, FS_LEGEND, LW = 26, 18, 16, 4
const PALETTE = [:orange, :seagreen, :dodgerblue, :firebrick, :royalblue]

_axis(fig_pos; kwargs...) = Axis(fig_pos; xlabelsize = FS_LABEL, ylabelsize = FS_LABEL,
    xticklabelsize = FS_TICK, yticklabelsize = FS_TICK, kwargs...)

# --- density overlay ---
fig = Figure(size = (900, 650))
ax = _axis(fig[1, 1]; xlabel = "r", ylabel = "ρ", yscale = log10)
for (i, (c, sol)) in enumerate(zip(cases, sols))
    lines!(ax, X, eval_density(sol, X); linewidth = LW, label = c.name, color = PALETTE[i])
end
axislegend(ax; labelsize = FS_LEGEND, position = :rt)
save("density_comparison.pdf", fig)

# --- energy breakdown, grouped by case ---
fields = (:Ekin, :Ecou, :Ehar, :Eexc, :Etot)
xs = Int[]
grp = Int[]
vals = Float64[]
for (fi, f) in enumerate(fields), (ci, sol) in enumerate(sols)
    push!(xs, fi)
    push!(grp, ci)
    push!(vals, getfield(sol.energies, f))
end

fig = Figure(size = (900, 650))
ax = _axis(fig[1, 1]; xlabel = "Energy component", ylabel = "Energy (Ha)",
    xticks = (1:length(fields), collect(string.(fields))))
barplot!(ax, xs, vals; dodge = grp, color = PALETTE[grp])
elems = [PolyElement(color = PALETTE[i]) for i in eachindex(cases)]
axislegend(ax, elems, [c.name for c in cases]; labelsize = FS_LEGEND, position = :rb)
save("energy_breakdown_comparison.pdf", fig)

# --- orbital energies, grouped by case ---
shells = sort(unique(occ[1] for sol in sols for occ in sol.occupied),
    by = s -> AtomicKohnSham.parse_shell(s)[1:2])
xs = Int[]
grp = Int[]
vals = Float64[]
for (si, shell) in enumerate(shells), (ci, sol) in enumerate(sols)
    m = findfirst(occ -> occ[1] == shell, sol.occupied)
    m === nothing && continue
    push!(xs, si)
    push!(grp, ci)
    push!(vals, sol.occupied[m][2])
end

fig = Figure(size = (900, 650))
ax = _axis(fig[1, 1]; xlabel = "Shell", ylabel = "Orbital energy (Ha)",
    xticks = (1:length(shells), shells))
barplot!(ax, xs, vals; dodge = grp, color = PALETTE[grp])
elems = [PolyElement(color = PALETTE[i]) for i in eachindex(cases)]
axislegend(ax, elems, [c.name for c in cases]; labelsize = FS_LEGEND, position = :rb)
save("orbital_energies_comparison.pdf", fig)
