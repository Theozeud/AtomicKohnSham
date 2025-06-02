include("../../../benchmarktools/atoms/setup.jl")
using KohnShamResolution

# LOG CONFIG
logconfig = LogConfig()

# Sodium
Sodium      = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(11, 11), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "Na", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Magnesium 
Magnesium   = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(12, 12), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "Be", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# ALumininum 
Alumininum  = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(13, 13), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "Al", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Silicium 
Silicium    = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(14, 14), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "Si", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Phosphore
Phosphore   = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(15, 15), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "P", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Soufre
Soufre      = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(16, 16), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "S", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Chlore
Chlore      = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(17, 17), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "Cl", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# Argon
Argon       = AtomProblem(;
                T               = Float64, 
                lh              = 1, 
                method          = CDA(0.3), 
                model           = ReducedHartreeFock(18, 18), 
                Rmax            = 80.0, 
                Nmesh           = 80,
                typemesh        = geometricmesh, 
                optsmesh        = (s=0.9,), 
                typebasis       = P1IntLegendreGenerator, 
                optsbasis       = (ordermax = 5,), 
                name            = "Ar", 
                scftol          = 1e-10,
                maxiter         = 100,
                hartree         = true,
                logconfig       = logconfig)

# RESOLUTION
@time "Sodium"      sol_sodium      = groundstate(Sodium)
@time "Magnesium"   sol_magnesium   = groundstate(Magnesium)
@time "Alumininum"  sol_aluminium   = groundstate(Alumininum)
@time "Silicium"    sol_silicium    = groundstate(Silicium)
@time "Phosphore"   sol_phosphore   = groundstate(Phosphore)
@time "Soufre"      sol_soufre      = groundstate(Soufre)
@time "Chlore"      sol_chlore      = groundstate(Chlore)
@time "Argon"       sol_argon       = groundstate(Argon)

# Print Results
original_stdout = stdout
output_file = open("results/atoms/rHF/recap global/second_row.txt", "w")
redirect_stdout(output_file)


println("Sodium")
orb_na = sol_sodium.orbitals_energy
println("1s : $(orb_na[1,1]) ; 2s : $(orb_na[1,2]); 2p : $(orb_na[2,1]); 3s : $(orb_na[1,3])")

println("Magnesium")
orb_mg = sol_magnesium.orbitals_energy
println("1s : $(orb_mg[1,1]) ; 2s : $(orb_mg[1,2]); 2p : $(orb_mg[2,1]); 3s : $(orb_mg[1,3])")

println("Alumininum")
orb_al = sol_aluminium.orbitals_energy
println("1s : $(orb_al[1,1]) ; 2s : $(orb_al[1,2]); 2p : $(orb_al[2,1]); 3s : $(orb_al[1,3]); 3p : $(orb_al[2,2])")

println("Silicium")
orb_si = sol_silicium.orbitals_energy
println("1s : $(orb_si[1,1]) ; 2s : $(orb_si[1,2]); 2p : $(orb_si[2,1]); 3s : $(orb_si[1,3]); 3p : $(orb_si[2,2])")

println("Phosphore")
orb_p = sol_phosphore.orbitals_energy
println("1s : $(orb_p[1,1]) ; 2s : $(orb_p[1,2]); 2p : $(orb_p[2,1]); 3s : $(orb_p[1,3]); 3p : $(orb_p[2,2])")

println("Soufre")
orb_s = sol_soufre.orbitals_energy
println("1s : $(orb_s[1,1]) ; 2s : $(orb_s[1,2]); 2p : $(orb_s[2,1]); 3s : $(orb_s[1,3]); 3p : $(orb_s[2,2])")

println("Chlore")
orb_cl = sol_chlore.orbitals_energy
println("1s : $(orb_cl[1,1]) ; 2s : $(orb_cl[1,2]); 2p : $(orb_cl[2,1]); 3s : $(orb_cl[1,3]); 3p : $(orb_cl[2,2])")

println("Argon")
orb_ar = sol_argon.orbitals_energy
println("1s : $(orb_ar[1,1]) ; 2s : $(orb_ar[1,2]); 2p : $(orb_ar[2,1]); 3s : $(orb_ar[1,3]); 3p : $(orb_ar[2,2])")


redirect_stdout(original_stdout)
close(output_file)
println("Finished")