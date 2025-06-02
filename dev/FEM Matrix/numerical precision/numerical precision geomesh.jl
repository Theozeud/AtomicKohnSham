using KohnShamResolution
using LinearAlgebra
using GenericLinearAlgebra
using DoubleFloats

N = 4
ordermax = 20
T = Float64
Rmin = T(1)
Rmax = T(50)
m = geometricmesh(Rmin, Rmax, N; T = T, s = 0.9)
Fbasis = P1IntLegendreGenerator(m, T; ordermax = ordermax)
@time "FA" FA       = sparse_stiffness_matrix(Fbasis)
@time "FM₀" FM₀     = sparse_mass_matrix(Fbasis)
@time "FM₋₁" FM₋₁   = sparse_weight_mass_matrix(Fbasis, -1)
@time "FM₋₂" FM₋₂   = sparse_weight_mass_matrix(Fbasis, -2)
@time "FF" FF       = weight_mass_3tensor(Fbasis, Monomial(-1))

T = Float64
Rmin = T(1)
Rmax = T(50)
m = geometricmesh(Rmin, Rmax, N; T = T, s = 0.9)
Dbasis = P1IntLegendreGenerator(m, T; ordermax = ordermax)
@time "DA" DA       = sparse_stiffness_matrix(Dbasis)
@time "DM₀" DM₀     = sparse_mass_matrix(Dbasis)
@time "DM₋₁" DM₋₁   = sparse_weight_mass_matrix(Dbasis, -1)
@time "DM₋₂" DM₋₂   = sparse_weight_mass_matrix(Dbasis, -2)
@time "DF" DF       = weight_mass_3tensor(Dbasis, Monomial(-1))


@show norm(DA - FA)
@show norm(DM₀ - FM₀)
@show norm(DM₋₁ - FM₋₁)
@show norm(DM₋₂ - FM₋₂)
@show norm(DF  - FF)