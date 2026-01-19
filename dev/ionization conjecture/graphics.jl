using AtomicKohnSham
import AtomicKohnSham: solve!, L_QUANTUM_LABELS
using GLMakie
using CairoMakie

# ===================== STEP 1 : PARAMETERS =====================#
# ATOM
Z = 118
N = 118

# FEM
Rmax        = 10000
Nmesh       = 70
ordermax    = 15
lmax        = 5

# ODA
maxiter     = 100
degen_tol   = 1e-3
scftol      = 1e-9

# PLOT
Rmax_plot   = 100
N_plot      = 10000
lplot       = 2

# ===================== STEP 2 : SOLVE =====================#

function run_solve(Z::Real, N::Real,
    Rmax::Real, Nmesh::Int, ordermax::Int, lmax::Int,
    maxiter::Int, degen_tol::Real, scftol::Real)
    model = KSEModel(; z = Z, N = N)
    mesh = expmesh(0,Rmax,Nmesh; s=1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
    dis = KSEDiscretization(lmax, basis, mesh, 1, 10)
    solver = KSESolver(model, dis, ODA(0.4);
            scftol = scftol,
            maxiter = maxiter,
            degen_tol = degen_tol,
            verbose = 0)
    solve!(solver)
    KSESolution(solver, "Z=$Z, N=$N"), solver
end

sol, solver = run_solve(Z, N, Rmax, Nmesh, ordermax, lmax, maxiter, degen_tol, scftol)

# ===================== STEP 3 : EVALUATION =====================#

function compute_potentials(solver, lmax::Int, Rmax_plot::Real, N_plot::Int)
    @unpack N,z = solver.model
    X = LinRange(0.001,Rmax_plot, N_plot)
    # Hartree
    @unpack cache, basis, Rmax = solver.discretization
    Vhart = evaluate(basis, cache.tmp_C, X)./ X .+ N/Rmax
    # Kinetic Potential
    Vl = zeros(length(X),lmax+1)
    for l ∈ 0:lmax
        @views Vl_l = Vl[:,l+1]
        @. Vl_l = l*(l+1)/(2*X^2)
    end
    # Nucleous-Electrons Potential
    Vnuc = @. - z/X
    # Total Potential
    V = zeros(length(X),lmax+1)
    for l ∈ 0:lmax
        @views V_l = V[:,l+1]
        @views Vl_l = Vl[:,l+1]
        @. V_l = Vhart + Vl_l + Vnuc
    end
    Vhart, Vnuc, Vl, V
end

Vhart, Vnuc, Vl, V = compute_potentials(solver, lmax, Rmax_plot, N_plot)

function compute_orbitals(sol, Rmax_plot::Real, N_plot::Int)
    X = LinRange(0,Rmax_plot, N_plot)
    u = zeros(length(X), length(sol.occupation_number))
    for (i,orb) ∈ enumerate(sol.occupation_number)
        idx = orb[1]
        u[:,i] = eigenvector(sol, idx, X)
        M = max(u[:,i]...)
        m = min(u[:,i]...)
        if abs(M) > abs(m)
            u[:,i] .*= sign(M)
        else
            u[:,i] .*= sign(m)
        end
    end
    u
end

u = compute_orbitals(sol, Rmax_plot, N_plot)

# ===================== STEP 4 : GRAPHICS =====================#

function make_figure(u::AbstractMatrix{<:Real}, V::AbstractMatrix{<:Real}, lplot::Int,
                     Rmax_plot::Real, N_plot::Int)
    X = LinRange(0,Rmax_plot, N_plot)
    fig = Figure(size = (1000, 800), fontsize = 30)

    Energies = sol.datas.orbitals_energy[:,lplot+1]
    min_E = min(Energies...)

    ax = Axis(fig[1, 1],
        xlabel = L"\text{r}",
        #yscale = log10,
        #xscale = log10,
        xticklabelsize = 40,
        yticklabelsize = 40,
        xlabelsize = 50,
        ylabelsize = 50,

    )
    ylims!(ax, min_E*1.2, 2)
    plots = []
    labels = []

    # Plot Potential
    V_l = V[:,lplot+1]
    li = lines!(ax, X, V_l; linewidth = 4)
    push!(plots, li)
    push!(labels, "V")

    # Plots Orbitals
    Colors = [:green, :orange, :purple, :yellow, :cyan]
    jc = 0
    for (i,orb) ∈ enumerate(sol.occupation_number)
        idx = orb[1]
        l = findfirst(x->x==idx[2],L_QUANTUM_LABELS) - 1
        if l == lplot
            jc += 1
            E = orb[2]
            hlines!(ax, E, linewidth=4, linestyle=:dash, label="E", color=Colors[jc])
            ui = u[:,i] .+ E
            li = lines!(ax, X, ui; linewidth = 4, color=Colors[jc])
            push!(plots, li)
            push!(labels, idx)
        end
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

X = LinRange(0.001,Rmax_plot, N_plot)

CairoMakie.activate!()
for l ∈ 0:lmax
    fig = make_figure(u, V, l, Rmax_plot, N_plot)
    save("./dev/ionization conjecture/fig_l=$l.pdf", fig)
end
