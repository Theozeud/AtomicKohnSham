using KohnShamResolution

include("../../benchmarktools/hydrogenoid/setup.jl")

using DoubleFloats


const proborder2 = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 90, 
                            Nmesh         = 10,
                            typemesh      = expmesh, 
                            typebasis     = P1IntLegendreGenerator, 
                            optsmesh      = (s = 2.0,),  
                            optsbasis     = (ordermax = 2, ),                           
                            name          = "IntLeg2-expmesh",
                            nU            = nothing)

const proborder5 = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 90, 
                            Nmesh         = 10,
                            typemesh      = expmesh, 
                            typebasis     = P1IntLegendreGenerator, 
                            optsmesh      = (s = 2.0,),  
                            optsbasis     = (ordermax = 5, ),                           
                            name          = "IntLeg5-expmesh",
                            nU            = nothing)

const proborder10 = HydrogenoidProblem(; 
                            T             = Float64, 
                            z             = 1, 
                            l             = 0, 
                            Rmax          = 90, 
                            Nmesh         = 5,
                            typemesh      = expmesh, 
                            typebasis     = P1IntLegendreGenerator, 
                            optsmesh      = (s = 1.0,),  
                            optsbasis     = (ordermax = 20, ),                           
                            name          = "IntLeg20-expmesh",
                            nU            = nothing)
#@time λgeo = eigvals_hydro(problinmesh)

#@time sol = eigen_hydro(problinmesh)
##plt = plot_eigenvector(1, Ugeo, probgeomesh)

convNmesh = convergenceNmesh(2 .^(3:6), [proborder2, proborder5, proborder10]; nums = [1])

convergence_plot_Nmesh(convNmesh)