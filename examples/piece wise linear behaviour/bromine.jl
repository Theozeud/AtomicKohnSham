using AtomicKohnSham
using Libxc


function ground_state_energy_carbon(N::AbstractVector)

    E_groundstate = zeros(Float64, length(N))

    for i ∈ eachindex(N)

        # Problem
        problem = AtomProblem(;
                                T = Float64,
                                z = 35,
                                N = N[i],
                                hartree = 1,
                                ex = Functional(:lda_x, n_spin = 2),
                                ec = Functional(:lda_c_pw, n_spin = 2),
                                lh = 2,
                                nh = 4,
                                Nmesh = 15,
                                Rmax = 400,
                                typemesh = expmesh,
                                optsmesh = (s = 1.5,),
                                typebasis = P1IntLegendreBasis,
                                optsbasis = (ordermax = 10,),
                                integration_method = GaussLegendre,
                                optsintegration = (npoints = 500,),
                                alg = ODA(0.4),
                                scftol = 1e-6,
                                maxiter = 100,
                                degen_tol = 1e-2,
                                logconfig = LogConfig(),
                                verbose = 0)

        # RESOLUTION
        @time sol = groundstate(problem);

        # SAVE ENERGY
        E_groundstate[i] = sol.energies[:Etot]

    end
    E_groundstate
end

using CairoMakie
using LaTeXStrings

N = LinRange(34,36,21)
E_groundstate = ground_state_energy_carbon(N)

fig = Figure(resolution = (1500, 1200))
ax = Axis(fig[1,1],
    xlabel = L"N",
    ylabel = L"E_{35,N}^{opt}",
    xticklabelsize = 55,
    yticklabelsize = 45,
    xlabelsize = 65,
    ylabelsize = 70
)
l1 = lines!(ax, N, E_groundstate, linewidth = 13, color = :royalblue)
l2 = begin
    coeff_34 = (E_groundstate[11] - E_groundstate[1])
    true_Egroundstate_34 = coeff_34 * (N[1:11] .- 34) .+ E_groundstate[1]

    coeff_35 = (E_groundstate[end] - E_groundstate[11])
    true_Egroundstate_35 = coeff_35 * (N[11:21] .- 35) .+ E_groundstate[11]

    true_Egroundstate = vcat(true_Egroundstate_34[1:end-1], true_Egroundstate_35)
    lines!(ax, N, true_Egroundstate, linestyle = (:dash,:dense), linewidth = 13, color = :tomato)
end

# Légende avec Legend
Legend(
        fig[1, 1],
        [l1, l2],
        ["Computed", "Expected"],
        halign = :right,
        valign = :top,
        tellheight = false,
        tellwidth = false,
        margin = (20, 20, 20, 20),
        orientation = :vertical,
        labelsize = 60,
        patchsize = (30, 20)
    )

save("./examples/piece wise linear behaviour/bromine.pdf",fig)
