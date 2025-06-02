include("../../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig()

# Lithium
Lithium     = AtomProblem(;
                T               = Float64, 
                lh              = 0, 
                method          = CDA(0.4), 
                model           = ReducedHartreeFock(3, 3), 
                Rmax            = 40.0, 
                Nmesh           = 60,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 4,), 
                name            = "Li", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Beryllium 
Beryllium   = AtomProblem(;
                T               = Float64, 
                lh              = 0, 
                method          = CDA(0.4), 
                model           = ReducedHartreeFock(4, 4), 
                Rmax            = 40.0, 
                Nmesh           = 60,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 4,), 
                name            = "Be", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Bore 
Bore        = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.4), 
                model           = ReducedHartreeFock(5, 5), 
                Rmax            = 40.0, 
                Nmesh           = 60,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 4,), 
                name            = "B", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Carbon 
Carbon      = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.4), 
                model           = ReducedHartreeFock(6, 6), 
                Rmax            = 40.0, 
                Nmesh           = 60,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "C", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Azote
Azote       = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.4), 
                model           = ReducedHartreeFock(7, 7), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "N", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Oxygen
Oxygen      = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(z, N), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "O", 
                scftol          = 1e-8,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Fer
Fer         = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(9, 9), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "F", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Neon
Neon        = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(10, 10), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "Ne", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# RESOLUTION
@time "Lithium"     sol_lithium     = groundstate(Lithium)
@time "Beryllium"   sol_beryllium   = groundstate(Beryllium)
@time "Bore"        sol_bore        = groundstate(Bore)
@time "Carbon"      sol_carbon      = groundstate(Carbon)
@time "Azote"       sol_azote       = groundstate(Azote)
@time "Oxygen"      sol_oxygen      = groundstate(Oxygen)
@time "Fer"         sol_fer         = groundstate(Fer)
@time "Neon"        sol_neon        = groundstate(Neon)

# Print Results
original_stdout = stdout
output_file = open("results/atoms/rHF/recap global/first_rows.txt", "w")
redirect_stdout(output_file)


println("Lithium")
orb_li = sol_lithium.orbitals_energy
println("1s : $(orb_li[1,1]) ; 2s : $(orb_li[1,2])")

println("Beryllium")
orb_be = sol_beryllium.orbitals_energy
println("1s : $(orb_be[1,1]) ; 2s : $(orb_be[1,2])")

println("Bore")
orb_b = sol_bore.orbitals_energy
println("1s : $(orb_b[1,1]) ; 2s : $(orb_b[1,2]) ; 2p : $(orb_b[2,1])")

println("Carbon")
orb_c = sol_carbon.orbitals_energy
println("1s : $(orb_c[1,1]) ; 2s : $(orb_c[1,2]) ; 2p : $(orb_c[2,1])")

println("Azote")
orb_n = sol_azote.orbitals_energy
println("1s : $(orb_n[1,1]) ; 2s : $(orb_n[1,2]) ; 2p : $(orb_n[2,1])")

println("Oxygen")
orb_o = sol_oxygen.orbitals_energy
println("1s : $(orb_o[1,1]) ; 2s : $(orb_o[1,2]) ; 2p : $(orb_o[2,1])")

println("Fer")
orb_f = sol_fer.orbitals_energy
println("1s : $(orb_f[1,1]) ; 2s : $(orb_f[1,2]) ; 2p : $(orb_f[2,1])")

println("Neon")
orb_ne = sol_neon.orbitals_energy
println("1s : $(orb_ne[1,1]) ; 2s : $(orb_ne[1,2]) ; 2p : $(orb_ne[2,1])")


redirect_stdout(original_stdout)
close(output_file)
println("Finished")