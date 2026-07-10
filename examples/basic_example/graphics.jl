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





# ===================================================================
#                           Plot Density
# ===================================================================
function plot_density(sol::KSESolution, X)
    fig = Figure(size = (1500, 1200))

    ax = Axis(fig[1, 1],
        xlabel = L"\text{r}",
        ylabel = L"\rho",
        #xscale = log10,
        yscale = log10,
        xticklabelsize = 50,
        yticklabelsize = 50,
        xlabelsize = 60,
        ylabelsize = 60
    )

    plots = []
    labels = []

    ρX = eval_density(sol, X)
    l = lines!(ax, X, ρX; linewidth = 8)
    push!(plots, l)
    push!(labels, sol.name)


    Legend(
        fig[1, 1],
        plots,
        labels,
        halign = :right,
        valign = :top,
        tellheight = false,
        tellwidth = false,
        margin = (20, 20, 20, 20),
        orientation = :vertical,
        labelsize = 50,
        patchsize = (20, 20)
    )

    return fig
end



# ===================================================================
#                           Plot Potentials
# ===================================================================
