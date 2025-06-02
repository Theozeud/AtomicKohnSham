using KohnShamResolution
using LinearAlgebra


N = 50

T = Float64
Rmin = T(0)
Rmax = T(50)
m = linmesh(Rmin, Rmax, N; T = T)
basis = ShortP1IntLegendreBasis(m, T; normalize = true, ordermax = 4)
deriv_basis = deriv(basis)
@time A   = mass_matrix(deriv_basis)
@time M₀  = mass_matrix(basis)
@time M₋₁ = weight_mass_matrix(basis, -1)
@time M₋₂ = weight_mass_matrix(basis, -2)

@show cond(A)
@show cond(M₀)
@show cond(M₋₁)
@show cond(M₋₂)

λ, U = eigen(0.5*A - M₋₁, M₀)



