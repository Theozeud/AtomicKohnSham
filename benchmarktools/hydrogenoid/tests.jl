include("new.jl")


prob = HydrogenoidProblem(;
    T = Float64,
    Z = 0.01,
    l = 0,
    Rmax = 1000000,
    Nmesh = 250,
    mesh_type = expmesh,
    basis_type = P1IntLegendreBasis,
    mesh_opts = (s = 1.0,),
    basis_opts = (ordermax = 5,),
    name = "IntLeg2-expmesh",
    eigs_idx = collect(1:10))

sol = eigen_hydro(prob)

#=
H,A, M0, basis = assemble_hydro_operators(prob)

using CairoMakie
X = LinRange(0,70,1000)
Y = evaluate(basis,U[:,4],X)


fun = theoretical_eigenvector(prob,[5])[1]
funY= fun.(X)

l = lines(X,Y)
lines!(X,funY)


import AtomicKohnSham: FunWeight
weight = FunWeight(fun;is_vectorized=false,is_inplace=false)
V=mass_vector(basis;weight=weight)
C = inv(M0)*V
YV = evaluate(basis,C,X)
lines!(X,YV)
l

A   = Symmetric(stiffness_matrix(basis))
h1norm = (U[:,4]-C)' * (A+M0) * (U[:,4]-C)
=#
