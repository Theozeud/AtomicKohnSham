using AtomicKohnSham

include("../../benchmarktools/hydrogenoid/setup.jl")

const problinmesh = HydrogenoidProblem(;
    T = Float64,
    z = 1,
    l = 0,
    Rmax = 300,
    Nmesh = 70,
    typemesh = linmesh,
    typebasis = P1IntLegendreBasis,
    optsmesh = (),
    optsbasis = (ordermax = 2,),
    name = "IntLeg2-linmesh",
    nU = nothing)

const probgeomesh = HydrogenoidProblem(;
    T = Float64,
    z = 1,
    l = 1,
    Rmax = 300,
    Nmesh = 70,
    typemesh = geometricmesh,
    typebasis = P1IntLegendreBasis,
    optsmesh = (s = 0.9,),
    optsbasis = (ordermax = 2,),
    name = "IntLeg2-geomesh",
    nU = nothing)

const probexpmesh = HydrogenoidProblem(;
    T = Float64,
    z = 1,
    l = 1,
    Rmax = 80,
    Nmesh = 70,
    typemesh = expmesh,
    typebasis = P1IntLegendreBasis,
    optsmesh = (s = 1.0,),
    optsbasis = (ordermax = 2,),
    name = "IntLeg2-expmesh",
    nU = nothing)

const probexpmesh2 = HydrogenoidProblem(;
    T = Float64,
    z = 1,
    l = 1,
    Rmax = 80,
    Nmesh = 70,
    typemesh = expmesh,
    typebasis = P1IntLegendreBasis,
    optsmesh = (s = 1.5,),
    optsbasis = (ordermax = 2,),
    name = "IntLeg2-expmesh2",
    nU = nothing)
#@time λgeo = eigvals_hydro(problinmesh)

#@time sol = eigen_hydro(problinmesh)
##plt = plot_eigenvector(1, Ugeo, probgeomesh)

convNmesh = convergenceNmesh(
    2 .^ (3:8), [problinmesh, probexpmesh, probexpmesh2]; nums = [2])

convergence_plot_Nmesh(convNmesh)
