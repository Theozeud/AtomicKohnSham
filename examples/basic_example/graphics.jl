using LaTeXStrings

# ===================================================================
#                 Plot Stopping Criteria through iterations
# ===================================================================
function plot_stopping_criteria(sol::KSESolution)
    fig = Figure(size = (1500, 1200), fontsize = 30)

    ax = Axis(fig[1, 1];
        xlabel = L"\text{Iteration}",
        ylabel = L"\text{SCF stopping criterion}",
        yscale = log10,
        xticklabelsize = 40,
        yticklabelsize = 40,
        xlabelsize = 50,
        ylabelsize = 50
    )

    sc = sol.logbook.stopping_criteria
    iters = eachindex(sc)

    lines!(ax, iters, sc; linewidth = 6)
    scatter!(ax, iters, sc; markersize = 30)

    return fig
end


# ===================================================================
#                           Plot Orbitals
# ===================================================================
