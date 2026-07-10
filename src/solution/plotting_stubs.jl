# ===================================================================
#                    PLOTTING (CairoMakie extension)
# ===================================================================
# The actual implementations live in ext/AtomicKohnShamCairoMakieExt.jl and
# are only compiled in when the user has `CairoMakie` loaded. These
# catch-all fallbacks give a clear error instead of a generic MethodError
# when that extension isn't active; they're automatically shadowed by the
# extension's more specific methods once CairoMakie is loaded.

function _requires_cairomakie(name::Symbol)
    error("`$(name)` requires CairoMakie to be loaded first, e.g. `using CairoMakie`.")
end

"""
    plot_density(sol::KSESolution, X)
    plot_density(sols::AbstractVector{<:KSESolution}, X)

Plot the radial electron density (log scale) on the radial points `X`, for
one solution or several overlaid. Requires `using CairoMakie`.
"""
plot_density(args...; kwargs...) = _requires_cairomakie(:plot_density)

"""
    plot_orbitals(sol::KSESolution, X; l = nothing)

Plot the occupied Kohn–Sham orbitals, each shifted to its orbital energy and
overlaid on the effective potential of its angular momentum channel
(see [`eval_effective_potential`](@ref)). Restrict to a single channel with
`l`. Requires `using CairoMakie`.
"""
plot_orbitals(args...; kwargs...) = _requires_cairomakie(:plot_orbitals)

"""
    plot_potentials(sol::KSESolution, X; l = 0, σ = 1)

Plot the nuclear, Hartree, exchange–correlation, and centrifugal
contributions to the effective potential of channel `l`, together with their
sum. Requires `using CairoMakie`.
"""
plot_potentials(args...; kwargs...) = _requires_cairomakie(:plot_potentials)

"""
    plot_convergence(sol::KSESolution)

Plot the SCF stopping criterion and the total-energy convergence
(`|Etot_i - Etot_final|`) against iteration, both on a log scale. Requires
`using CairoMakie`.
"""
plot_convergence(args...; kwargs...) = _requires_cairomakie(:plot_convergence)

"""
    plot_energy_breakdown(sol::KSESolution)

Bar chart of the energy decomposition (`Ekin`, `Ecou`, `Ehar`, `Eexc`,
`Etot`). Requires `using CairoMakie`.
"""
plot_energy_breakdown(args...; kwargs...) = _requires_cairomakie(:plot_energy_breakdown)
