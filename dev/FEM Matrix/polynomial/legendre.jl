using KohnShamResolution
using Test
using Plots


T = Float64

for n ∈ 1:1
    P = intLegendre(n)
    println(P)
end

