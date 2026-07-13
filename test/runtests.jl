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
include("hydrogen.jl")
include("scandium.jl")
