using KohnShamResolution
using LinearAlgebra
using BenchmarkTools

using KohnShamResolution: fill_stiffness_matrix!, fill_mass_matrix! 

original_stdout = stdout
output_file = open("tests/Performance/fem/matrices_$(Threads.nthreads())_threads.txt", "a")
redirect_stdout(output_file)

println("")
println("=====================================")
println("=====================================")

# IntLeg2 - Nmesh = 50
println("IntLeg2, Nmesh = 50")
Nmesh = 50
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, Nmesh; T = T)
println("basis")
@time basis = P1IntLegendreGenerator(m, T; ordermax = 2)
l = length(basis)

println("A")
A = zeros(T,l,l)
@time fill_stiffness_matrix!(basis, A)
println("M0")
M₀ = zeros(T,l,l)
@time fill_mass_matrix!(basis, M₀)
println("M₋₁")
M₋₁ = zeros(T,l,l)
@time fill_mass_matrix!(basis, -1, M₋₁)
println("M₋₂")
M₋₂ = zeros(T,l,l)
@time fill_mass_matrix!(basis, -2, M₋₂)

println("=====================================")

# IntLeg5 - Nmesh = 50
println("IntLeg5, Nmesh = 50")
Nmesh = 50
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, Nmesh; T = T)
println("basis")
@time basis = P1IntLegendreGenerator(m, T; ordermax = 5)
l = length(basis)

println("A")
A = zeros(T,l,l)
@time fill_stiffness_matrix!(basis, A)
println("M0")
M₀ = zeros(T,l,l)
@time fill_mass_matrix!(basis, M₀)
println("M₋₁")
M₋₁ = zeros(T,l,l)
@time fill_mass_matrix!(basis, -1, M₋₁)
println("M₋₂")
M₋₂ = zeros(T,l,l)
@time fill_mass_matrix!(basis, -2, M₋₂)

println("=====================================")

# IntLeg10 - Nmesh = 300
println("IntLeg10, Nmesh = 300")
Nmesh = 300
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, Nmesh; T = T)
println("basis")
@time basis = P1IntLegendreGenerator(m, T; ordermax = 10)
l = length(basis)

println("A")
A = zeros(T,l,l)
@time fill_stiffness_matrix!(basis, A)
println("M0")
M₀ = zeros(T,l,l)
@time fill_mass_matrix!(basis, M₀)
println("M₋₁")
M₋₁ = zeros(T,l,l)
@time fill_mass_matrix!(basis, -1, M₋₁)
println("M₋₂")
M₋₂ = zeros(T,l,l)
@time fill_mass_matrix!(basis, -2, M₋₂)


redirect_stdout(original_stdout)
close(output_file)

println("Performance finished")