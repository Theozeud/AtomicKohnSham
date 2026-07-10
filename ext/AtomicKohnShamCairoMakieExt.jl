module AtomicKohnShamCairoMakieExt

using AtomicKohnSham
using CairoMakie

import AtomicKohnSham: parse_shell, has_exchcorr

# ===================================================================
#                       PUBLICATION STYLE
# ===================================================================
# Sizes tuned to stay readable once a figure is shrunk to journal column
# width; FIGSIZE is generous so labels/ticks/legend never crowd the plot.
const FS_TICK   = 34
const FS_LABEL  = 52
const FS_TITLE  = 40
const FS_LEGEND = 30
const LW        = 5
const FIGSIZE   = (1400, 1050)

# Colorblind-friendly palette, reordered so the first series drawn in any
# plot is never the plain, "default-looking" blue.
const _PALETTE = [:orange, :seagreen, :orchid, :dodgerblue, :firebrick, :royalblue, :gold]
_color(i::Int) = _PALETTE[mod1(i, length(_PALETTE))]

_orbital_l(shell::String) = parse_shell(shell)[2]

function _pub_axis(fig_pos; kwargs...)
    Axis(fig_pos; xlabelsize = FS_LABEL, ylabelsize = FS_LABEL,
        xticklabelsize = FS_TICK, yticklabelsize = FS_TICK,
        titlesize = FS_TITLE, kwargs...)
end

_pub_legend!(ax; position = :rt) =
    axislegend(ax; labelsize = FS_LEGEND, patchsize = (50, 20), position = position)

"""
Reasonable energy window for a given angular momentum channel: tight around
its occupied orbital energies rather than the full range of the effective
potential (which diverges to -∞ at r=0 and, for l>0, blows up there too).
Falls back to all occupied energies if the channel itself has none occupied.
"""
function _energy_ylims(sol::KSESolution, lval::Int)
    εs = [occ[2] for occ in sol.occupied if _orbital_l(occ[1]) == lval]
    isempty(εs) && (εs = [occ[2] for occ in sol.occupied])
    isempty(εs) && return nothing
    elo, ehi = extrema(εs)
    pad = max(0.3 * (ehi - elo), 0.15 * abs(elo), 0.2)
    return (elo - pad, max(ehi + pad, 0.2))
end

# ===================================================================
#                           DENSITY
# ===================================================================
function AtomicKohnSham.plot_density(sol::KSESolution, X::AbstractVector{<:Real})
    return AtomicKohnSham.plot_density([sol], X)
end

function AtomicKohnSham.plot_density(sols::AbstractVector{<:KSESolution},
                                    X::AbstractVector{<:Real})
    fig = Figure(size = FIGSIZE)
    ax = _pub_axis(fig[1, 1]; xlabel = "r", ylabel = "ρ", yscale = log10)

    for (i, sol) in enumerate(sols)
        ρX = eval_density(sol, X)
        lines!(ax, X, ρX; linewidth = LW, label = sol.name, color = _color(i))
    end

    _pub_legend!(ax; position = :rt)
    return fig
end

# ===================================================================
#                           ORBITALS
# ===================================================================
function AtomicKohnSham.plot_orbitals(sol::KSESolution, X::AbstractVector{<:Real};
                                     l::Union{Nothing, Int} = nothing)
    ls = l === nothing ? sort(unique(_orbital_l(occ[1]) for occ in sol.occupied)) : [l]

    fig = Figure(size = (FIGSIZE[1], FIGSIZE[2] * length(ls)))

    for (row, lval) in enumerate(ls)
        ax = _pub_axis(fig[row, 1]; xlabel = "r", ylabel = "Energy (Ha)",
            title = "l = $(lval)", xscale = log10)

        Veff = eval_effective_potential(sol, lval, X)
        lines!(ax, X, Veff; linewidth = LW, color = :black, label = "Effective potential")

        ylims = _energy_ylims(sol, lval)
        yspan = ylims === nothing ? 1.0 : (ylims[2] - ylims[1])

        jc = 0
        for occ in sol.occupied
            shell, ϵ, _ = occ
            _orbital_l(shell) == lval || continue
            jc += 1
            color = _color(jc)
            u = eval_orbital(sol, shell, X)
            # Orbitals are drawn for their shape, not their true (FEM-normalized)
            # amplitude: rescale each to a fixed fraction of the panel height so
            # tightly-bound (large |ε|) and loosely-bound orbitals are equally
            # visible on the same energy axis.
            umax = maximum(abs, u)
            scale = umax > 0 ? (0.08 * yspan) / umax : one(umax)
            hlines!(ax, [ϵ]; linewidth = 2, linestyle = :dash, color = color)
            lines!(ax, X, u .* scale .+ ϵ; linewidth = LW, color = color, label = shell)
        end

        ylims === nothing || ylims!(ax, ylims...)

        _pub_legend!(ax; position = :rb)
    end

    return fig
end

# ===================================================================
#                           POTENTIALS
# ===================================================================
function AtomicKohnSham.plot_potentials(sol::KSESolution, X::AbstractVector{<:Real};
                                       l::Int = 0, σ::Int = 1)
    model = sol.context.model

    fig = Figure(size = FIGSIZE)
    ax = _pub_axis(fig[1, 1]; xlabel = "r", ylabel = "Potential (Ha)", xscale = log10)

    i = 1
    lines!(ax, X, eval_nuclear(sol, X); linewidth = LW, label = "Nuclear", color = _color(i))
    i += 1
    lines!(ax, X, eval_kinetic_potential(sol, l, X); linewidth = LW,
        label = "Centrifugal (l=$(l))", color = _color(i))
    i += 1

    if !iszero(model.hartree)
        lines!(ax, X, eval_hartree(sol, X); linewidth = LW, label = "Hartree", color = _color(i))
        i += 1
    end

    if has_exchcorr(model)
        Vxc = model.nspin == 1 ? eval_vxc(sol, X) : eval_vxc(sol, X, σ)
        xclabel = model.nspin == 1 ? "Exchange-correlation" : "Exchange-correlation (σ=$(σ))"
        lines!(ax, X, Vxc; linewidth = LW, label = xclabel, color = _color(i))
        i += 1
    end

    lines!(ax, X, eval_effective_potential(sol, l, X; σ = σ); linewidth = LW + 1,
        label = "Total (effective)", color = :black)

    ylims = _energy_ylims(sol, l)
    ylims === nothing || ylims!(ax, ylims...)

    _pub_legend!(ax; position = :rb)
    return fig
end

# ===================================================================
#                           CONVERGENCE
# ===================================================================
function AtomicKohnSham.plot_convergence(sol::KSESolution{T}) where T
    fig = Figure(size = FIGSIZE)
    ax = _pub_axis(fig[1, 1]; xlabel = "Iteration", ylabel = "Convergence", yscale = log10)

    sc = sol.logbook.stopping_criteria
    iters = eachindex(sc)
    lines!(ax, iters, sc; linewidth = LW, label = "Stopping criterion", color = _color(1))

    Etot_final = sol.energies.Etot
    dE = [max(abs(E - Etot_final), eps(T)) for E in sol.logbook.Etot]
    lines!(ax, eachindex(dE), dE; linewidth = LW, label = "|Etot - Etot(final)|",
        color = _color(2))

    _pub_legend!(ax; position = :rt)
    return fig
end

# ===================================================================
#                       ENERGY BREAKDOWN
# ===================================================================
function AtomicKohnSham.plot_energy_breakdown(sol::KSESolution)
    labels = ["Ekin", "Ecou", "Ehar", "Eexc", "Etot"]
    values = [sol.energies.Ekin, sol.energies.Ecou, sol.energies.Ehar,
        sol.energies.Eexc, sol.energies.Etot]

    fig = Figure(size = FIGSIZE)
    ax = _pub_axis(fig[1, 1]; xticks = (1:5, labels), ylabel = "Energy (Ha)")
    barplot!(ax, 1:5, values; color = [_color(i) for i in 1:5])

    return fig
end

end # module
