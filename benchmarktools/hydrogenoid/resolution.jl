#############################################################################################
# RESOLUTION
#############################################################################################

# Normamization procedure
function L2normalization(U, M₀)
    ax = axes(U, 1)
    nU = zero(U)
    for k in axes(U, 2)
        normalization = sqrt(sum([U[i, k] * U[j, k] * M₀[i, j] for i in ax for j in ax]) *
                             4π)
        nU[:, k] .= U[:, k] / normalization
    end
    nU
end

# Eigenvalues and Eigenvectors
function eigen_hydro(problem)
    @unpack T, z, l, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis = problem
    mesh = typemesh(zero(T), Rmax, Nmesh; T = T, optsmesh...)
    basis = typebasis(mesh, T; optsbasis...)
    A = Symmetric(stiffness_matrix(basis))
    M₀ = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(mass_matrix(basis, -1))
    if l == 0
        H = T(0.5) * A - z*M₋₁
        Λ, U = eigen(H, M₀)
        nU = L2normalization(U, M₀)
        return HydrogenoidSolution(problem, Λ, nU)
    else
        M₋₂ = Symmetric(mass_matrix(basis, -2))
        H = T(0.5) * (A + l * (l+1) * M₋₂) - z * M₋₁
        Λ, U = eigen(H, M₀)
        nU = L2normalization(U, M₀)
        return HydrogenoidSolution(problem, Λ, nU)
    end
end

# Eigenvalue with Eigvals (usefull for Double64)
function eigvals_hydro(problem)
    @unpack T, z, l, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis = problem
    mesh = typemesh(zero(T), Rmax, Nmesh; T = T, optsmesh...)
    basis = typebasis(mesh, T; optsbasis...)
    A = Symmetric(stiffness_matrix(basis))
    M₀ = Symmetric(mass_matrix(basis))
    M₋₁ = Symmetric(mass_matrix(basis, -1))
    if l == 0
        H = T(0.5) * A - z * M₋₁
        Λ = eigvals(inv(M₀) * H)
        return HydrogenoidSolution(problem, Λ, [1])
    else
        M₋₂ = Symmetric(mass_matrix(basis, -2))
        H = T(0.5) * (A + l * (l+1) * M₋₂) - z * M₋₁
        Λ = eigvals(inv(M₀) * H)
        return HydrogenoidSolution(problem, Λ, [1])
    end
end
