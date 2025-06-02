using KohnShamResolution

original_stdout = stdout
output_file = open("tests/Performance/basis.txt", "w")
redirect_stdout(output_file)

# Création des éléments puis de la base
T = Float64
mesh = linmesh(0,10,1000)

# Tests Creation of elements
println("Elements")
@time intLeg = IntLegendreElements(T; ordermin = 2, ordermax = 2, binf = -T(1), bsup = T(1))

# Tests Creation of the basis
println("IntLegBasis")
@time basintleg = ShortIntLegendreBasis(mesh, T; ordermin = 2, ordermax = 5)

# Tests Creation of NewBasis
println("IntLegGenerator")
@time new_basis = IntLegendreGenerator(mesh, T; ordermin = 2, ordermax = 5)

redirect_stdout(original_stdout)
close(output_file)
println("Performance finished")