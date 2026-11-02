cd(@__DIR__)
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using AtomicKohnSham
using CairoMakie
using Printf

# Abaque: |Etot - Etot_ref| vs Rmax (log-log), for Z=50 (tin, Z=N=50),
# nspin=1, Slater exchange ONLY -- ec = NoFunctional(1), no correlation. One
# curve per (Nmesh, ordermax) combination: Nmesh via color (a single-hue
# sequential ramp, since Nmesh is an ordered quantity), ordermax via line
# style.
#
# Same study and grid as "examples/abaque/oxygen" / "examples/abaque/Z20",
# generalized to Z=50. Two things differ, both because tin is considerably
# harder to converge than the lighter atoms:
#
# - Rmax_ref = 1500, not 2000: at Rmax=2000 the SCF stalls at MAXITERS for
#   tin even with a generous (Nmesh, ordermax) (checked directly: 300
#   iterations, still not converged), while Rmax=1500 with the same
#   discretization converges cleanly in under 100 iterations. Not fully
#   understood why 2000 specifically is harder -- flagged here rather than
#   chased further, since 1500 already sits comfortably beyond the grid's
#   own Rmax=1000 ceiling.
# - target_err = 1e-2, much looser than the lighter atoms: like calcium (see
#   "examples/abaque/Z20"), the grid's own Rmax=1000 ceiling carries real
#   truncation error for tin -- here roughly 3.5e-3 Ha, an order of
#   magnitude worse than calcium's, consistent with heavier atoms needing
#   bigger simulation domains. Not a reference-quality issue -- see STEP 1.

# ===================== STEP 1 : PARAMETERS =====================#
Z, N = 50, 50
basis_size(Z) = (lh = min(ceil(Int, Z / 20), 3), nh = clamp(Z, 2, 10))
lh, nh = basis_size(Z)

ex = BuiltinFunctional(:lda_x; nspin = 1)
ec = NoFunctional(1)
model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)

Nmesh_values    = [16, 32,64]
ordermax_values = [6, 10, 14]
Rmax_values     = [100, 200, 500, 1000, 2000, 5000]

scftol  = 1e-8
maxiter = 200

# Reference: comfortably beyond every grid corner in all three knobs at
# once. Rmax_ref=1500, not 2000 -- see note above STEP 1.
Rmax_ref     = 5000
Nmesh_ref    = 70
ordermax_ref = 14
scftol_ref   = 1e-8
maxiter_ref  = 400

target_err = 1e-2   # grid's own Rmax=1000 ceiling carries ~3.5e-3 Ha of
                     # truncation error for tin -- see note above STEP 1.

# ===================== STEP 2 : SOLVER HELPER =====================#
"""
Solve the Z=N=50, nspin=1, Slater-exchange-only (no correlation) ground state
at the given discretization. Returns `(success, Etot)`, or `nothing` if the
run errors out entirely.
"""
function run_case(Rmax::Real, Nmesh::Int, ordermax::Int, scftol::Real, maxiter::Int)
    try
        mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
        basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
        dis = KSEDiscretization(basis, model; lh = lh, nh = nh)

        aufbau = OptimizedAufbau(max_degen = 2)
        alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = scftol)
        sol = groundstate(model, dis, alg; maxiter = maxiter)

        return (success = sol.success, Etot = sol.energies.Etot)
    catch e
        @warn "run_case failed" Rmax Nmesh ordermax exception = e
        return nothing
    end
end

# ===================== STEP 3 : REFERENCE =====================#
ref = run_case(Rmax_ref, Nmesh_ref, ordermax_ref, scftol_ref, maxiter_ref)
ref === nothing && error("reference calculation failed -- cannot build the abaque")
Etot_ref = ref.Etot
println("Reference (", ref.success, "): Etot = ", Etot_ref,
    "  (Rmax=", Rmax_ref, ", Nmesh=", Nmesh_ref, ", ordermax=", ordermax_ref, ")")

# ===================== STEP 4 : GRID =====================#
# err[i, k, j] = error at (Nmesh_values[i], ordermax_values[k], Rmax_values[j])
err     = fill(NaN, length(Nmesh_values), length(ordermax_values), length(Rmax_values))
outcome = fill("", length(Nmesh_values), length(ordermax_values), length(Rmax_values))

for (i, Nmesh) in enumerate(Nmesh_values), (k, ordermax) in enumerate(ordermax_values),
        (j, Rmax) in enumerate(Rmax_values)
    c = run_case(Rmax, Nmesh, ordermax, scftol, maxiter)
    if c === nothing
        outcome[i, k, j] = "FAILED"
        continue
    end
    err[i, k, j] = abs(c.Etot - Etot_ref)
    outcome[i, k, j] = c.success
end

finite_err = filter(!isnan, err)
best_err = isempty(finite_err) ? NaN : minimum(finite_err)
println("Best grid error vs reference: ", best_err)
if isnan(best_err) || best_err > target_err
    println("WARNING: best grid point only reaches ", best_err,
        " -- did not get below ", target_err,
        ". Either the reference (Rmax_ref=", Rmax_ref, ", Nmesh_ref=", Nmesh_ref,
        ", ordermax_ref=", ordermax_ref, ") needs pushing further, OR the grid's own",
        " Rmax=", maximum(Rmax_values), " ceiling is genuinely not enough to reach",
        " target_err for this atom -- check which by comparing against a larger Rmax_ref.")
else
    println("OK: best grid point reaches ", target_err, " -- reference is well resolved.")
end

# ===================== STEP 5 : ABAQUE.TXT =====================#
open("abaque.txt", "w") do io
    println(io, "Z=50 (tin, Z=N=50), nspin=1, Slater exchange only (no correlation)")
    println(io, "|Etot - Etot_ref| vs (Nmesh, Rmax), one block per ordermax")
    println(io, "="^90)
    println(io)
    println(io, "Reference: Etot = ", Etot_ref, "  (Rmax=", Rmax_ref, ", Nmesh=", Nmesh_ref,
        ", ordermax=", ordermax_ref, ", scftol=", scftol_ref, ")")
    println(io, "Best grid error vs reference: ", best_err,
        best_err <= target_err ? "  (<= target $target_err)" : "  (ABOVE target $target_err)")
    println(io)
    colw = 13
    for (k, ordermax) in enumerate(ordermax_values)
        println(io, "ordermax = ", ordermax)
        println(io, "-"^90)
        println(io, rpad("Nmesh \\ Rmax", 14), join(rpad(string(R), colw) for R in Rmax_values))
        for (i, Nmesh) in enumerate(Nmesh_values)
            row = join(rpad(isnan(err[i, k, j]) ? "FAILED" : @sprintf("%.4e", err[i, k, j]), colw)
                        for j in eachindex(Rmax_values))
            println(io, rpad(string(Nmesh), 14), row)
        end
        println(io)
    end
    println(io, "SUCCESS/FAILURE per point")
    println(io, "="^90)
    for (k, ordermax) in enumerate(ordermax_values)
        println(io, "ordermax = ", ordermax)
        println(io, "-"^90)
        println(io, rpad("Nmesh \\ Rmax", 14), join(rpad(string(R), colw) for R in Rmax_values))
        for (i, Nmesh) in enumerate(Nmesh_values)
            row = join(rpad(outcome[i, k, j], colw) for j in eachindex(Rmax_values))
            println(io, rpad(string(Nmesh), 14), row)
        end
        println(io)
    end
end

# ===================== STEP 6 : PLOT =====================#
const FS_LABEL, FS_TICK, FS_LEGEND = 24, 16, 14
_axis(fig_pos; kwargs...) = Axis(fig_pos; xlabelsize = FS_LABEL, ylabelsize = FS_LABEL,
    xticklabelsize = FS_TICK, yticklabelsize = FS_TICK, kwargs...)
_floor_err(v) = max(v, 1e-17)

# Color = Nmesh (ordered quantity -> single-hue sequential ramp, light to
# dark, steps 250/500/700 of the blue ramp in the dataviz skill's reference
# palette). Line style = ordermax (a second, independent knob).
const NMESH_COLORS    = ["#86b6ef", "#256abf", "#0d366b"]
const ORDERMAX_STYLES = [:solid, :dash, :dot]

fig = Figure(size = (1150, 700))
ax = _axis(fig[1, 1]; xlabel = "Rmax", ylabel = "|Etot - Etot_ref| (Ha)",
    xscale = log10, yscale = log10,
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9),
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(9),
    title = "Z=50, nspin=1, Slater exchange only", titlesize = FS_LABEL)
for (i, Nmesh) in enumerate(Nmesh_values), (k, ordermax) in enumerate(ordermax_values)
    ys = [_floor_err(err[i, k, j]) for j in eachindex(Rmax_values)]
    lines!(ax, Rmax_values, ys; linewidth = 3, linestyle = ORDERMAX_STYLES[k], color = NMESH_COLORS[i])
    scatter!(ax, Rmax_values, ys; markersize = 8, color = NMESH_COLORS[i])
end

color_elems = [LineElement(color = c, linewidth = 3) for c in NMESH_COLORS]
color_labels = ["Nmesh=$n" for n in Nmesh_values]
style_elems = [LineElement(color = :black, linestyle = s, linewidth = 3) for s in ORDERMAX_STYLES]
style_labels = ["ordermax=$o" for o in ordermax_values]
Legend(fig[1, 2], color_elems, color_labels, "Nmesh"; labelsize = FS_LEGEND, titlesize = FS_LEGEND)
Legend(fig[1, 3], style_elems, style_labels, "ordermax"; labelsize = FS_LEGEND, titlesize = FS_LEGEND)

save("abaque_Z50.pdf", fig)
