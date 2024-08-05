using Plots

function compute_sol(;z, N, Rmax, Nmesh, lₕ, maxiter, oda, tol = 1e-5, T = Float64)
    method = ConstantODA(T(oda))
    m = linmesh(zero(T), Rmax, Nmesh)
    KM = KohnShamExtended(z = z, N = N)

    basis_p1 = ShortP1Basis(m, T; left = false, right = false, normalize = true)
    basis_intleg2 = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = 2, normalize = true)
    basis_intleg3 = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = 3, normalize = true)
    basis_intleg4 = ShortP1IntLegendreBasis(m, T; left = false, right = false, ordermin = 2, ordermax = 4, normalize = true)

    # Final Discretization
    D1 = KohnShamRadialDiscretization(lₕ, basis_p1, m)
    D2 = KohnShamRadialDiscretization(lₕ, basis_intleg2, m)
    D3 = KohnShamRadialDiscretization(lₕ, basis_intleg3, m)
    D4 = KohnShamRadialDiscretization(lₕ, basis_intleg4, m)

    # Solution
    @time "With P1" sol1 = groundstate(KM, D1, method; tol = tol, hartree = true, maxiter = maxiter, potential = :pde)
    @time "With IntLeg2" sol2 = groundstate(KM, D2, method; tol = tol, hartree = true, maxiter = maxiter, potential = :pde)
    @time "With IntLeg3" sol3 = groundstate(KM, D3, method; tol = tol, hartree = true, maxiter = maxiter, potential = :pde)
    @time "With IntLeg4" sol4 = groundstate(KM, D4, method; tol = tol, hartree = true, maxiter = maxiter, potential = :pde)
    sols = [sol1, sol2, sol3, sol4]
    # Title
    title = "Rmax = $Rmax, z = $z, oda = $oda, N = $N, Nmesh = $Nmesh"
    # Label
    label = ["P1","IntLeg2","IntLeg3", "IntLeg4"]
    # Return
    sols, title, label
end


function plot_crit(sols, tol; title, label)
    plt_criteria = plot( size = (650,500), margin = 0.5Plots.cm, legend = :outertopright, yaxis=:log,
    legendfontsize  = 12,  
    titlefontsize   = 12,
    guidefontsize   = 12,
    tickfontsize    = 12)
    xlabel!(plt_criteria, "Iteration")
    ylabel!(plt_criteria, "Convergence criteria")
    title!(plt_criteria, title)
    lmax = 0
    for (i,sol) ∈ enumerate(sols)
        plot!(plt_criteria, sol.crit, lw = 4, label = label[i], markershape = :x, markersize = 10)
        if length(sol.crit)>lmax
            lmax = length(sol.crit)
        end
        end
    plot!(ones(lmax)*tol, label = "tolerance", ls = :dash, color = :black, lw =4 )
    plt_criteria
end

function plot_Ehisto(sols; title, label)
    plt_Ehisto = plot(  size = (650,500), margin = 0.5Plots.cm, legend = :outertopright,
                        legendfontsize  = 11,  
                        titlefontsize   = 11,
                        guidefontsize   = 11,
                        tickfontsize    = 11)
    xlabel!(plt_Ehisto, "Iteration")
    ylabel!(plt_Ehisto, "Evolution of ϵ₁")
    title!(plt_Ehisto, title)
    for (i,sol) ∈ enumerate(sols)
        plot!(plt_Ehisto, sol.Ehisto, lw = 3, label = label[i], markershape = :x, markersize = 10)
    end
    plt_diffE = plot(  size = (650,500), margin = 0.5Plots.cm, legend = :outertopright, yaxis = :log,
                        legendfontsize  = 12,  
                        titlefontsize   = 12,
                        guidefontsize   = 12,
                        tickfontsize    = 12)
    xlabel!(plt_diffE, "Iteration")
    ylabel!(plt_diffE, "Convergence of ϵ₁")
    title!(plt_diffE, title)
    for (i,sol) ∈ enumerate(sols)
        plot!(plt_diffE, eachindex(sol.crit)[2:end], abs.(sol.Ehisto[2:end] .- sol.Ehisto[1:end-1]), lw = 3, label = label[i], markershape = :x, markersize = 10)
    end
    plt_Ehisto, plt_diffE
end

function plot_density(sols, Rmax; title, label)
    plt_ρ = plot(  size = (650,500), margin = 0.5Plots.cm, legend = :outertopright,
                        legendfontsize  = 11,  
                        titlefontsize   = 11,
                        guidefontsize   = 11,
                        tickfontsize    = 11)
    xlabel!(plt_ρ, "r")
    ylabel!(plt_ρ, "ρ")
    title!(plt_ρ, title)
    X = range(0, Rmax, 1000)
    for (i,sol) ∈ enumerate(sols)
        plot!(plt_ρ, X, sol.ρ.(X), lw = 3, label = label[i])
    end
    plt_ρ
end