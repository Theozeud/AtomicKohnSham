using KohnShamResolution
using Plots
using LinearAlgebra

using GenericLinearAlgebra
using DoubleFloats

# Theoretical Eigenvalue and Eigenvectors
vect_theo(n, Rmax, T) = x-> sqrt(T(2)/T(Rmax)) * sin(n*T(π)*x/T(Rmax)) 
val_theo(n, Rmax, T)  = (T(π)/T(Rmax))^2 * n^2

##############################################################################
# Function to compute eigenvalue and eigenvector for different discretizations

# With finite difference
function eigen_with_finite_diff(Nmesh, Rmax, T)
    A = SymTridiagonal(2*ones(T, Nmesh), - ones(T, Nmesh-1)) .* ((T(Nmesh) - T(1))/ T(Rmax))^2
    eigvals(A)
end

# With P1
function eigen_with_P1(mesh, T)
    normalize = true
    left = false
    right = false
    basis = ShortP1Basis(mesh, T; left = left, right = right, normalize = normalize)
    M₀  = mass_matrix(basis)
    AP1   = Symmetric(mass_matrix(deriv(basis)))
    eigvals(inv(M₀) * AP1)
end

# With P1-Integrated Legendre Polynomials
function eigen_with_P1IntLeg(mesh, ordermax, T)
    normalize = true
    ordermin = 2
    left = false
    right = false
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = ordermin, ordermax = ordermax,  normalize = normalize, left = left, right = right)
    M₀  = mass_matrix(basis)
    AIntleg2   = Symmetric(mass_matrix(deriv(basis)))
    eigvals(inv(M₀) * AIntleg2)
end

# With P1
function eigen_with_P1(mesh, n, T)
    normalize = true
    left = false
    right = false
    Rmax = T(last(mesh))
    basis = ShortP1Basis(mesh, T; left = left, right = right, normalize = normalize)
    M₀  = Symmetric(mass_matrix(basis))
    AP1   = Symmetric(mass_matrix(deriv(basis)))
    eign = build_on_basis(basis, eigen(AP1, M₀).vectors[:,n])
    eign / sqrt(integrate(eign*eign, T(Rmin), T(Rmax)))
end

# With P1-Integrated Legendre Polynomials
function eigen_with_P1IntLeg(mesh, ordermax, n, T)
    normalize = true
    ordermin = 2
    left = false
    right = false
    Rmax = T(last(mesh))
    basis = ShortP1IntLegendreBasis(mesh, T; ordermin = ordermin, ordermax = ordermax,  normalize = normalize, left = left, right = right)
    M₀  = Symmetric(mass_matrix(basis))
    AIntleg2   = Symmetric(mass_matrix(deriv(basis)))
    eign = build_on_basis(basis, eigen(AIntleg2, M₀).vectors[:,n])
    eign/ sqrt(integrate(eign*eign, T(Rmin), T(Rmax)))
end

##############################################################################
# Plot Eigenvalue

function plot_eigenvalue(Nmesh, Rmax, order = [2,3], T)
    # Creation of the mesh
    m = linmesh(Rmin, Rmax, Nmesh, T)
    
    # With finite difference
    #vals_diff = eigen_with_finite_diff(Nmesh, Rmax, T)
    # With P1
    vals_p1 = eigen_with_P1(m, T)
    # With P1-Integrated Legendre Polynomials
    vals_il = []
    for o ∈ order
        vals_ilo = eigen_with_P1IntLeg(m, o, T)
        push!(vals_il, vals_ilo)
    end

    plt_eigenvalue_intleg = plot(size= (1500, 700), margin = 0.5Plots.cm,                
                                 legendfontsize  = 14,  
                                 titlefontsize   = 14,
                                 guidefontsize   = 14,
                                 tickfontsize    = 14)
    plot!(vals_p1, label = "P1", lw = 4)
    max_index = length(vals_p1)
    for (i,o) ∈ enumerate(order)
        plot!(vals_il[i], label = "P1 + Integrated Legendre ordre "*string(o), lw = 4)
        max_index = max(length(vals_p1), length(vals_il[i]))
    end
    plot!(val_theo.(1:max_index, Rmax), label = "Theoretical", lw = 4)
    return plt_eigenvalue_intleg
end

##############################################################################
# Plot Eigenvector

function plot_eigenvector(Nmesh, Rmax, n, order)
    # Creation of the mesh
    m = linmesh(Rmin, Rmax, Nmesh, T)
   
    # With finite difference
    #_, eigs_diff = eigen_with_finite_diff(Nmesh, Rmax)
    # plot!(plt_eigenvector, range(Rmin, Rmax, Nmesh), eigs_diff[:,n]/ sqrt((Rmax - Rmin)/(Nmesh-1)), label = "Finite Difference", lw = 2)
    # With P1
    #eigs_p1 = eigen_with_P1(m, n)
    # plot!(plt_eigenvector, X, (sign(vect_theo(n, Rmax)(1) * eigs_p1(1)) * eigs_p1).(X), label = "P1", lw = 2)
       
    # With Integrated Legendre 
    eigs_il = eigen_with_P1IntLeg(m, order, n, T)
    
    X = LinRange(Rmin, Rmax, 10000)
    plt_eigenvector = plot( size= (1500, 700), margin = 1Plots.cm,                
                            legendfontsize  = 14,  
                            titlefontsize   = 14,
                            guidefontsize   = 14,
                            tickfontsize    = 14)
    xlabel!("r")
    title!(string(n)*"-th eigenvector")
    plot!(plt_eigenvector, X, vect_theo(n, Rmax).(X), label = "theo n = "*string(n), lw = 4)
    plot!(plt_eigenvector, X, (sign(vect_theo(n, Rmax)(1) * eigs_il(1)) * eigs_il).(X), label = "P1 + Integrated Legendre ordre "*string(order), lw = 4)
    
    plt_eigenvector
end

##################################################################################
# Compute Error

function compute_error_eigenvalue(vecNmesh, Rmax, num, T)

    label = ["Différence finie", "P1", "IntLeg 2", "IntLeg 3", "IntLeg 4"]
    ϵerror = zeros(T, length(vecNmesh), length(label))
    for (i, Nmesh) ∈ enumerate(vecNmesh)
        # Creation of the mesh
        m = linmesh(Rmin, Rmax, Nmesh, T)
        
        # True eigenvalue
        @time "Nmesh = "*string(Nmesh) begin
        true_val = val_theo.(1:Nmesh, Rmax, T)

        # With finite difference
        vals_diff = eigen_with_finite_diff(Nmesh, Rmax)

        # With P1
        vals_p1 = eigen_with_P1(m, T)

        # With P1-Integrated Legendre Polynomials ordre 2
        vals_il2 = eigen_with_P1IntLeg(m, 2, T)

        # With P1-Integrated Legendre Polynomials ordre 3
        vals_il3 = eigen_with_P1IntLeg(m, 3, T)

        # With P1-Integrated Legendre Polynomials ordre 4
        vals_il4 = eigen_with_P1IntLeg(m, 4, T)
        end

        # Compute the error for eigenvalues and the fundamental
        ϵerror[i,1] = abs(vals_diff[num] - true_val[num])
        ϵerror[i,2] = abs(vals_p1[num] - true_val[num])
        ϵerror[i,3] = abs(vals_il2[num] - true_val[num])
        ϵerror[i,4] = abs(vals_il3[num] - true_val[num])
        ϵerror[i,5] = abs(vals_il4[num] - true_val[num])
    end

    # Creation of the plot for eigenvalue
    plt_ϵ_error = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :topright, xaxis=:log, yaxis=:log,
                        legendfontsize  = 14,  
                        titlefontsize   = 14,
                        guidefontsize   = 14,
                        tickfontsize    = 14)
    xlabel!(plt_ϵ_error, "Nmesh")
    ylabel!(plt_ϵ_error, "Error on the $num-th eigenvalue")
    title!(plt_ϵ_error, "Rmax = $Rmax")
    for i ∈ eachindex(label)
        plot!(plt_ϵ_error, vecNmesh, ϵerror[:, i], lw = 4, label = label[i], markershape = :x, markersize = 10)
    end
    plt_ϵ_error
end

##############################################

#=
# General Discretization Parameters
T = Double64
Rmin = zero(T)
vecNmesh = 2 .^(3:7)
=#

#=
plt_eigenvalue = plot_eigenvalue(100, 10, [2,3])

savefig(plt_eigenvalue, "test/image_tests/Laplacien - Valeurs Propres 1-2-3")
plt_eigenvalue = plot_eigenvalue(100,10, [2,3,4,5])
savefig(plt_eigenvalue, "test/image_tests/Laplacien - Valeurs Propres 1-2-3-4-5")


plt_eigenvector = plot_eigenvector(100,10, 198, 3)
savefig(plt_eigenvector, "test/image_tests/Laplacien - Vecteur Propre 198 - ordre 3")
plt_eigenvector = plot_eigenvector(100,10, 199, 3)
savefig(plt_eigenvector, "test/image_tests/Laplacien - Vecteur Propre 199 - ordre 3")
plt_eigenvector = plot_eigenvector(100,10, 200, 3)
savefig(plt_eigenvector, "test/image_tests/Laplacien - Vecteur Propre 200 - ordre 3")

plt_eigenvalue_error = compute_error_eigenvalue(vecNmesh, 10, 1)
savefig(plt_eigenvalue_error, "test/image_tests/Laplacien - Erreurs valeurs propres")

compute_error_eigenvalue(vecNmesh, 10, 1)
=#
