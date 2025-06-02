using KohnShamResolution
using Plots
using LinearAlgebra

# Utils to compute second membre
using Integrals
function compute_second_membre(f, basis)
    T = KohnShamResolution.bottom_type(basis)
    b = zeros(T, length(basis))
    for i ∈ eachindex(basis)
        (ib, iib) = KohnShamResolution.find_basis(basis, i)
        spb = KohnShamResolution.getbasis(basis, ib)
        val = zero(T)
        for (j,_,_, invϕ) ∈ spb.infos[iib]
            Pelement = KohnShamResolution.getnormalization(spb,iib) * getpolynomial(spb.elements, j)
            g(x) = f(invϕ(x)) * Pelement(x)
            prob = IntegralProblem((x,p) -> g(x), -T(1), T(1))
            val += invϕ[1] * solve(prob, QuadGKJL(); abstol = 1e-13).u
        end
        b[i] = val
    end
    b
end

# With P1
function solve_p1(mesh, f, T, power = -1)
    left = false
    right = false
    normalize = true
    basis = ShortP1Basis(mesh, T; left = left, right = right, normalize = normalize)
    CP1   = weight_mass_matrix(basis, power)
    #F = [integrate(f * build_basis(basis, i), Rmin, Rmax) for i ∈ eachindex(basis)]
    F = compute_second_membre(f, basis)
    build_on_basis(basis, CP1\F)
end

# With P1-Integrated Legendre Polynomials
function solve_intleg(mesh, ordermax, f, T, power = -1)
    left = false
    right = false
    ordermin = 2
    normalize = true
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = ordermin, ordermax = ordermax,  normalize = normalize, left = left, right = right)
    CIL   = weight_mass_matrix(basis, power)
    #F = [integrate(f * build_basis(basis, i), Rmin, Rmax) for i ∈ eachindex(basis)]
    F = compute_second_membre(f, basis)
    build_on_basis(basis, CIL\F)
end

#=
using LinearSolve
using IncompleteLU
using SparseArrays
using Preconditioners
function solve_intleg_ls(mesh, ordermax, f)
    left = false
    right = false
    ordermin = 2
    normalize = true
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = ordermin, ordermax = ordermax,  normalize = normalize, left = left, right = right)
    CIL   = weight_mass_matrix(basis, -1)
    F = [integrate(f * build_basis(basis, i), Rmin, Rmax) for i ∈ eachindex(basis)]
    prob = LinearProblem(CIL, F)
    sol = solve(prob)
    build_on_basis(basis, sol)
end
=#


# Plot solutions
function plot_sol(mesh, f, T, power = -1)
    f_true = Monomial(-power, T(1)) * f
    X = LinRange(Rmin, Rmax, 1000)
    plt = plot(size = (1000, 800), margin = 0.5Plots.cm)
    plot!(plt, X, solve_p1(mesh, f, T, power).(X),  label = "p1", lw = 3)
    plot!(plt, X, solve_intleg(mesh, 2, f, T, power).(X), label = "Integrated Legendre ordre 2", lw = 3)
    plot!(plt, X, solve_intleg(mesh, 3, f, T, power).(X), label = "Integrated Legendre ordre 3", lw = 3)
    plot!(plt, X, f_true.(X), label = "Theoretical", lw = 2, color = :black, ls = :dash)
    title!("Nmesh = $(length(mesh))")
    xlabel!("r")
    plt
end

# Compute error
function plot_error(vecNmesh, f, Rmin, Rmax, T = Float64, power = -1)
    f_true = Monomial(-power, T(1)) * f
    label = ["P1", "IntLeg 2", "IntLeg 3", "IntLeg4", "IntLeg5"]
    ϵerror = zeros(T, length(vecNmesh), length(label))
    for (i, Nmesh) ∈ enumerate(vecNmesh)

        # Creation of the mesh
        m = linmesh(Rmin, Rmax, Nmesh; T = T)
        
        # With P1
        sol_p1 = solve_p1(m, f, T, power)

        # With P1-Integrated Legendre Polynomials ordre 2
        sol_il2 = solve_intleg(m, 2, f, T, power)

        # With P1-Integrated Legendre Polynomials ordre 3
        sol_il3 = solve_intleg(m, 3, f, T, power)

        # With P1-Integrated Legendre Polynomials ordre 4
        sol_il4 = solve_intleg(m, 4, f, T, power)

        # With P1-Integrated Legendre Polynomials ordre 5
        sol_il5 = solve_intleg(m, 5, f, T, power)

        # Compute the error for eigenvalues and the fundamental
        c = (T(Rmax)-T(Rmin))/(Nmesh - 1)
        ϵerror[i,1] = sqrt(c * sum(abs.(sol_p1.(m.points[1:end-1])  .- f_true.(m.points[1:end-1])).^2))
        ϵerror[i,2] = sqrt(c * sum(abs.(sol_il2.(m.points[1:end-1]) .- f_true.(m.points[1:end-1])).^2))
        ϵerror[i,3] = sqrt(c * sum(abs.(sol_il3.(m.points[1:end-1]) .- f_true.(m.points[1:end-1])).^2))
        ϵerror[i,4] = sqrt(c * sum(abs.(sol_il4.(m.points[1:end-1]) .- f_true.(m.points[1:end-1])).^2))
        ϵerror[i,5] = sqrt(c * sum(abs.(sol_il5.(m.points[1:end-1]) .- f_true.(m.points[1:end-1])).^2))
    end

    # Creation of the plot for eigenvalue
    plt_error = plot(size = (1000,800), margin = 0.5Plots.cm, legend = :topright, xaxis=:log, yaxis=:log,
                     legendfontsize  = 14,  
                     titlefontsize   = 14,
                     guidefontsize   = 14,
                     tickfontsize    = 14)
    xlabel!(plt_error, "Nmesh")
    ylabel!(plt_error, "L2 Error")
    title!(plt_error, "Rmax = $Rmax")
    for i ∈ eachindex(label)[1:5]
        plot!(plt_error, vecNmesh, ϵerror[:, i], lw = 4, label = label[i], markershape = :x, markersize = 10)
    end
    plt_error
end


#=
# General Discretization Parameters
T = Float64
Rmin = 0.5
Rmax = 1
Nmesh = 300
m = linmesh(Rmin, Rmax,Nmesh)
=#

#=
X = LinRange(0, 0.001, 10000)
m = linmesh(Rmin, Rmax, 10)
plot(X, solve_intleg(m, 4).(X), label = "Integrated Legendre ordre 4", lw = 3, size = (1000,500))
m = linmesh(Rmin, Rmax, 100)
plot!(X, solve_intleg(m, 4).(X), label = "Integrated Legendre ordre 4", lw = 3, size = (1000,500))
plot!(X, f_true.(X), label = "Theoretical", lw = 2, color = :black, ls = :dash)
=#

#=
f = -Polynomial([Rmin*Rmax, -Rmin - Rmax, 1], 0)
plt_error = plot_error(2 .^(2:7), f)
savefig(plt_error, "test/image_tests/Coulomb - Error degf2 1-2-3-4-5")

f = -Polynomial([Rmin*Rmax, -Rmin - Rmax, 1], 3)
plt_error = plot_error(2 .^(2:7), f)
savefig(plt_error, "test/image_tests/Coulomb - Error degf5 1-2-3-4-5")
=#