using AtomicKohnSham
using Libxc
using Test
using LinearAlgebra
using SparseArrays
using JLD2

#=
using Aqua
Aqua.test_all(AtomicKohnSham)
=#

include("fem.jl")
include("utils.jl")
include("discretization.jl")
include("physics.jl")
include("algorithms.jl")
include("solver.jl")
include("solution.jl")
include("hydrogen.jl")
include("scandium.jl")
