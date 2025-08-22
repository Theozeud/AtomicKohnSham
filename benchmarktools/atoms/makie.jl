function plot_stopping_criteria(sols)
    fig = Figure(size = (1500, 1200), fontsize = 30)

    ax = Axis(fig[1, 1],
        xlabel = L"\text{Iteration}",
        ylabel = L"\text{SCF tolerance}",
        yscale = log10,
        xticklabelsize = 40,
        yticklabelsize = 40,
        xlabelsize = 50,
        ylabelsize = 50
    )

    plots = []
    labels = []

    for sol in sols
        sc = sol.log.stopping_criteria
        iters = 1:length(sc)
        push!(plots, lines!(ax, iters, sc; linewidth = 6))
        push!(plots, scatter!(ax, iters, sc; markersize = 30))
        push!(labels, sol.name)
    end

    legend = [[plots[2 * i - 1], plots[2 * i]]=>(;) for i in eachindex(sols)]
    Legend(
        fig[1, 1],
        legend,
        labels,
        halign = :right,
        valign = :top,
        tellheight = false,
        tellwidth = false,
        margin = (20, 20, 20, 20),
        orientation = :vertical
    )

    return fig
end

function plot_density(sols, X)
    fig = Figure(size = (1500, 1200), fontsize = 30)

    ax = Axis(fig[1, 1],
        xlabel = L"\text{r}",
        ylabel = L"\rho",
        #xscale = log10,
        yscale = log10,
        xticklabelsize = 40,
        yticklabelsize = 40,
        xlabelsize = 50,
        ylabelsize = 50
    )

    plots = []
    labels = []

    for sol in sols
        mesh = sol.discretization.mesh.points
        ρX = eval_density(sol, X)
        l = lines!(ax, X, ρX; linewidth = 6)
        push!(plots, l)
        push!(labels, sol.name)
        #scatter!(ax, mesh, eval_density(sol, mesh); markersize = 30)

    end

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

function plot_eigenvector(sol, X)
    fig = Figure(size = (1500, 1200), fontsize = 30)

    ax = Axis(fig[1, 1],
        xlabel = L"\text{r}",
        ylabel = L"\epsilon",
        yscale = log10,
        #scale = log10,
        xticklabelsize = 40,
        yticklabelsize = 40,
        xlabelsize = 50,
        ylabelsize = 50
    )

    plots = []
    labels = []

    mesh = sol.discretization.mesh.points
    for orb ∈ sol.occupation_number
        idx = orb[1]
        ϕnlσ = abs.(eigenvector(sol, idx, X)) ./ X
        l = lines!(ax, X, ϕnlσ; linewidth = 6)
        push!(plots, l)
        push!(labels, idx)
        scatter!(ax, mesh, abs.(eigenvector(sol, idx, mesh)) ./ mesh; markersize = 30)
    end


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
