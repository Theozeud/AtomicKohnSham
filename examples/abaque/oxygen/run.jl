cd(@__DIR__)
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using AtomicKohnSham
using CairoMakie
using Printf

# Abaque: |Etot - Etot_ref| vs Rmax (log-log), for oxygen (Z=N=8), nspin=1,
# Slater exchange ONLY -- ec = NoFunctional(1), no correlation. One curve per
# (Nmesh, ordermax) combination: Nmesh via color (a single-hue sequential
# ramp, since Nmesh is an ordered quantity), ordermax via line style.
#
# Same study as "examples/abaque/hydrogen", generalized to Z=8. Two things
# differ from the hydrogen version, both driven by oxygen actually having
# more electrons to place:
#
# - Nmesh_values starts at 8, not 5: `basis_size(8) = (lh=1, nh=8)` needs at
#   least 8 radial basis functions per l-channel just to HAVE 8 orbitals to
#   fill -- Nmesh=5/ordermax=2 (hydrogen's smallest grid point, dof=7) throws
#   `ArgumentError("illegal choice of eigenvalue indices...")` for any atom
#   needing nh > 7. Nmesh=8/ordermax=2 (dof=13) clears that floor.
# - ordermax_values stops at 10, same ceiling as hydrogen and for the same
#   reason (see hydrogen's run.jl for the ordermax=20+ instability details).
#
# Etot_ref comes from a single reference run, one clear notch beyond every
# grid corner in all three knobs at once. If the best (finest) grid corner
# doesn't get down near `target_err`, that is itself the signal the
# reference needs to be pushed further, not that the FEM discretization has
# a floor there -- see the check in STEP 4.

# ===================== STEP 1 : PARAMETERS =====================#
Z, N = 8, 8
basis_size(Z) = (lh = min(ceil(Int, Z / 20), 3), nh = clamp(Z, 2, 10))
lh, nh = basis_size(Z)

ex = BuiltinFunctional(:lda_x; nspin = 1)
ec = NoFunctional(1)
model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)

Nmesh_values    = [8, 16, 32]
ordermax_values = [2, 5, 10]
Rmax_values     = [2, 5, 10, 20, 50, 100, 300, 1000]

scftol  = 1e-8
maxiter = 200

# Reference: comfortably beyond every grid corner in all three knobs at once.
Rmax_ref     = 2000
Nmesh_ref    = 48
ordermax_ref = 13
scftol_ref   = 1e-9
maxiter_ref  = 300

target_err = 1e-8   # looser than hydrogen's 1e-10 -- an 8-electron SCF just
                     # doesn't settle as tightly as a 1-electron one; grid-best
                     # reaches ~9e-9 here (see description.md)

# ===================== STEP 2 : SOLVER HELPER =====================#
"""
Solve the Z=N=8, nspin=1, Slater-exchange-only (no correlation) ground state
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
    println(io, "Oxygen (Z=N=8), nspin=1, Slater exchange only (no correlation)")
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
    title = "O, nspin=1, Slater exchange only", titlesize = FS_LABEL)
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

save("abaque_oxygen.pdf", fig)
