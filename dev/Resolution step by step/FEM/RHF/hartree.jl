using KohnShamResolution
using Plots
using LinearAlgebra
using Integrals

######################################################################################################
# Function to compute the potential

# Formula with Integrals
function approx_intH10(mesh, ρ, T)
    Rmin = T(first(mesh))
    Rmax = T(last(mesh))
    f(x) =  x * ρ(x)
    g(x) = x^2 * ρ(x)
    vect = zeros(T, length(mesh))
    vect[1] =  KohnShamResolution.approximate_integral(f, (Rmin, Rmax); method = QuadGKJL(), reltol = 100*eps(T), abstol = 100*eps(T)) *Rmin
    for i ∈ eachindex(vect)[2:end]
        domain1 = (Rmin, T(mesh[i]))
        vect[i] = KohnShamResolution.approximate_integral(g, domain1; method = QuadGKJL(), reltol = 100*eps(T), abstol = 100*eps(T))
        domain2 = (T(mesh[i]), Rmax)
        vect[i] +=  mesh[i] * KohnShamResolution.approximate_integral(f, domain2; method = QuadGKJL(), reltol = 100*eps(T), abstol = 100*eps(T))
    end
    4π * (vect .- mesh ./ Rmax .* last(vect))
end

function approx_int(mesh, ρ, T)
    Rmin = first(mesh)
    Rmax = last(mesh)
    f(x) =  x * ρ(x)
    g(x) = x^2 * ρ(x)
    vect = zeros(T, length(mesh))
    vect[1] =  KohnShamResolution.approximate_integral(f, (Rmin, Rmax); method = QuadGKJL(), reltol = 100*eps(T), abstol = 100*eps(T))
    for i ∈ eachindex(vect)[2:end]
        domain1 = (first(mesh), mesh[i])
        vect[i] = 1/mesh[i] * KohnShamResolution.approximate_integral(g, domain1; method = QuadGKJL(), reltol = 100*eps(T), abstol = 100*eps(T))
        domain2 = (mesh[i], Rmax)
        vect[i] +=  KohnShamResolution.approximate_integral(f, domain2; method = QuadGKJL(), reltol = 100*eps(T), abstol = 100*eps(T))
    end
    4π * vect
end

# Resolution by PDE
# Computation of the potential given a basis
function potentialH10(basis, ρ)
    deriv_basis = deriv(basis)
    CP1   = mass_matrix(deriv_basis)
    f(x) = 4π*ρ(x) * x#Monomial(1)
    F =  weight_mass_vector(basis, f)
    poly = build_on_basis(basis, CP1\F)
    poly
end

function potential(basis, ρ, Rmin, Rmax)
    sol = potentialH10(basis, ρ)
    g = ρ * Monomial(2)
    Cᵨ = 4π * integrate(g, Rmin, Rmax)
    sol = sol * Monomial(-1) + Cᵨ / (Rmax-Rmin)
    sol
end

function potential_funH10(basis, ρ)
    deriv_basis = deriv(basis)
    CP1   = mass_matrix(deriv_basis)
    f(x) = 4π*ρ(x) * x
    F =  weight_mass_vector(basis, f)
    build_on_basis(basis, CP1\F)
end

function potential_fun(basis, ρ, Rmin, Rmax)
    sol = potentialH10(basis, ρ)
    g = ρ(x) * x^2
    Cᵨ = 4π * integrate(g, Rmin, Rmax)
    sol = sol * Monomial(-1) + Cᵨ / (Rmax-Rmin)
    sol
end

# With P1
function potential_p1H10(mesh, ρ, T)
    basis = ShortP1Basis(mesh, T; left = false, right = false, normalize = true)
    potentialH10(basis, ρ)
end

function potential_p1(mesh, ρ, T, Rmin, Rmax)
    basis = ShortP1Basis(mesh, T; left = false, right = false, normalize = true)
    potential(basis, ρ, Rmin, Rmax)
end

function potential_fun_p1H10(mesh, ρ, T)
    basis = ShortP1Basis(mesh, T; left = false, right = false, normalize = true)
    potential_funH10(basis, ρ)
end

function potential_fun_p1(mesh, ρ, T, Rmin, Rmax)
    basis = ShortP1Basis(mesh, T; left = false, right = false, normalize = true)
    potential_fun(basis, ρ, Rmin, Rmax)
end

# With P1-Integrated Legendre Polynomials
function potential_intlegH10(mesh, ordermax, ρ, T)
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = 2, ordermax = ordermax,  normalize = true, left = false, right = false)
    potentialH10(basis, ρ)
end

function potential_intleg(mesh, ordermax, ρ, T, Rmin, Rmax)
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = 2, ordermax = ordermax,  normalize = true, left = false, right = false)
    potential(basis, ρ, Rmin, Rmax)
end

function potential_fun_intlegH10(mesh, ordermax, ρ, T)
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = 2, ordermax = ordermax,  normalize = true, left = false, right = false)
    potential_funH10(basis, ρ)
end

function potential_fun_intleg(mesh, ordermax, ρ, T, Rmin, Rmax)
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = 2, ordermax = ordermax,  normalize = true, left = false, right = false)
    potential_fun(basis, ρ, Rmin, Rmax)
end

######################################################################################################
# Compute error for potential
function plot_errorH10(vecNmesh, ρ, Rmin, Rmax, T = Float64)
    label = ["P1", "IntLeg 2", "IntLeg 3", "IntLeg4", "IntLeg5"]
    ϵerror = zeros(T, length(vecNmesh), length(label))
    X = range(Rmin, Rmax, 1000)
    for (i, Nmesh) ∈ enumerate(vecNmesh)
        # Creation of the mesh
        m = linmesh(Rmin, Rmax, Nmesh; T = T)
        # Exact solution
        f_true = approx_intH10(X, ρ, T)
        # With P1
        sol_p1 = potential_p1H10(m, ρ, T)
        # With P1-Integrated Legendre Polynomials ordre 2
        sol_il2 = potential_intlegH10(m, 2, ρ, T)
        # With P1-Integrated Legendre Polynomials ordre 3
        sol_il3 = potential_intlegH10(m, 3, ρ, T)
        # With P1-Integrated Legendre Polynomials ordre 4
        sol_il4 = potential_intlegH10(m, 4, ρ, T)
        # With P1-Integrated Legendre Polynomials ordre 5
        sol_il5 = potential_intlegH10(m, 5, ρ, T)
        # Compute the error for eigenvalues and the fundamental
        c = (T(Rmax)-T(Rmin))/(length(X) - 1)
        ϵerror[i,1] = sqrt(c * sum(abs.(sol_p1.(X[2:end])  .- f_true[2:end]).^2))
        ϵerror[i,2] = sqrt(c * sum(abs.(sol_il2.(X[2:end]) .- f_true[2:end]).^2))
        ϵerror[i,3] = sqrt(c * sum(abs.(sol_il3.(X[2:end]) .- f_true[2:end]).^2))
        ϵerror[i,4] = sqrt(c * sum(abs.(sol_il4.(X[2:end]) .- f_true[2:end]).^2))
        ϵerror[i,5] = sqrt(c * sum(abs.(sol_il5.(X[2:end]) .- f_true[2:end]).^2))
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
    ϵerror, plt_error
end

function plot_error(vecNmesh, ρ, Rmin, Rmax, T = Float64)
    label = ["P1", "IntLeg 2", "IntLeg 3", "IntLeg4", "IntLeg5"]
    ϵerror = zeros(T, length(vecNmesh), length(label))
    X = range(Rmin, Rmax, 1000)
    for (i, Nmesh) ∈ enumerate(vecNmesh)
        # Creation of the mesh
        m = linmesh(Rmin, Rmax, Nmesh; T = T)
        # Exact solution
        f_true = approx_int(X, ρ, T)
        # With P1
        sol_p1 = potential_p1(m, ρ, T, Rmin, Rmax)
        # With P1-Integrated Legendre Polynomials ordre 2
        sol_il2 = potential_intleg(m, 2, ρ, T, Rmin, Rmax)
        # With P1-Integrated Legendre Polynomials ordre 3
        sol_il3 = potential_intleg(m, 3, ρ, T, Rmin, Rmax)
        # With P1-Integrated Legendre Polynomials ordre 4
        sol_il4 = potential_intleg(m, 4, ρ, T, Rmin, Rmax)
        # With P1-Integrated Legendre Polynomials ordre 5
        sol_il5 = potential_intleg(m, 5, ρ, T, Rmin, Rmax)

        # Compute the error for eigenvalues and the fundamental
        c = (T(Rmax)-T(Rmin))/(length(X) - 1)
        ϵerror[i,1] = sqrt(c * sum(abs.(sol_p1.(X[2:end])  .- f_true[2:end]).^2))
        ϵerror[i,2] = sqrt(c * sum(abs.(sol_il2.(X[2:end]) .- f_true[2:end]).^2))
        ϵerror[i,3] = sqrt(c * sum(abs.(sol_il3.(X[2:end]) .- f_true[2:end]).^2))
        ϵerror[i,4] = sqrt(c * sum(abs.(sol_il4.(X[2:end]) .- f_true[2:end]).^2))
        ϵerror[i,5] = sqrt(c * sum(abs.(sol_il5.(X[2:end]) .- f_true[2:end]).^2))
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
    ϵerror, plt_error
end



######################################################################################################
# Function to compute ther Hartree matrix
function hartree_v1(basis, ρ, Rmin, Rmax)
    potent(x) = potential(basis, ρ, Rmin, Rmax)(x)
    weight_mass_matrix(basis, potent)
end

function hartree_v2(basis, ρ, Rmin, Rmax)
    deriv_basis = deriv(basis)
    CIL   = mass_matrix(deriv_basis)
    f(x) = 4π*ρ(x) * x
    @show F =  weight_mass_vector(basis, f)
    @show coeff = CIL\F
    g = ρ * Monomial(2)
    @show Cᵨ = 4π * integrate(g, Rmin, Rmax)
    @show vectorweight_mass_matrix(basis, coeff, Monomial(-1))
    vectorweight_mass_matrix(basis, coeff, Monomial(-1))+ Cᵨ/(Rmax-Rmin) * mass_matrix(basis), F
end

function hartree_v3(basis, D, Rmin, Rmax)
    deriv_basis = deriv(basis)
    M₀ = mass_matrix(basis)
    A = mass_matrix(deriv_basis)
    F = weight_mass_3tensor(basis, Monomial(1))
    @show @tensor B[m] := D[i,j] * F[i,j,m]
    @show Pot = A\B
    @tensor MV[i,j] := Pot[k] * F[i,j,k]
    mw = weight_mass_matrix(basis, 2)
    @show @tensor Cᵨ = D[i,j] * mw[i,j]
    @show MV
    MV + Cᵨ/(Rmax-Rmin) * M₀, B
end


