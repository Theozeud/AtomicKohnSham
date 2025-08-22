# Convergence Plot with Nmesh

function convergence_plot_Nmesh(sols::HydrogenoidConvergenceNmesh; nums = (first(sols.num),))
    fig = Figure(resolution = (1300, 1000))
    ax = Axis(fig[1, 1];
        title = "Convergence Plot",
        xlabel = "Nmesh",
        ylabel = "Error",
        xscale = log10,
        yscale = log10,
        titlesize = 35,
        xlabelsize = 25,
        ylabelsize = 25,
        xticklabelsize = 30,
        yticklabelsize = 30
    )

    for prob in sols.probs
        for num in nums
            pos = findfirst(==(num), nums)
            scatterlines!(ax, sols.vecNmesh, sols.ΔΛ[prob.name][:, pos];
                linewidth = 4,
                label = "$(prob.name)-$num",
                marker = :x,
                markersize = 30
            )
        end
    end

    axislegend(ax, position = :lb, labelsize = 30)
    return fig
end


function convergence_plot_Rmax(sols::HydrogenoidConvergenceRmax; nums = (first(sols.num),))
    fig = Figure(resolution = (1300, 1000))
    ax = Axis(fig[1, 1];
        title = "Convergence Plot",
        xlabel = "Rmax",
        ylabel = "Error",
        xscale = log10,
        yscale = log10,
        titlesize = 12,
        xlabelsize = 12,
        ylabelsize = 12,
        xticklabelsize = 12,
        yticklabelsize = 12
    )

    for prob in sols.probs
        for num in nums
            scatterlines!(ax, sols.vecRmax, sols.ΔΛ[prob.name][:, num];
                linewidth = 4,
                label = "$(prob.name)-$num",
                marker = :x,
                markersize = 10
            )
        end
    end

    axislegend(ax, position = :lb, labelsize = 12)
    return fig
end
