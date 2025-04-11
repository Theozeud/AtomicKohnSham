using KohnShamResolution
using UnPack

using KohnShamResolution: NoWeight, default_method, find_intersection_indices!, getgenerator, getshift, getinvshift, getmesh, getprod, swsp,
AbstractWeight, IntegrationMethod, IntegrationData, InvX


function fill_mass_matrix2!( pb::PolynomialBasis, 
                            A::AbstractMatrix{<:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, precomputations = pb
    @show J = zeros(Int, 2, pb.max_length_intersection[1])
    for I ∈ pb.matrix_fill_indices[end-1:end]
        @show count = find_intersection_indices!(J,pb.indices_cells[I[1]],pb.indices_cells[I[2]])
        for c ∈ 1:count
            @show i = J[1,c]
            @show j = J[2,c]
            @show P = getgenerator(pb, I[1], i)
            @show Q = getgenerator(pb, I[2], j)
            @show ϕ = getshift(pb, I[1], i)
            @show invϕ = getinvshift(pb, I[1], i)
            @show (a,b) = getmesh(pb, I[1], i)
            iP = pb.indices_generators[I[1]][i]
            iQ = pb.indices_generators[I[2]][j]
            prod = getprod(pb, iP, iQ)
            @show intdata = IntegrationData(weight,
                        (P,Q),
                        prod,
                        ϕ,
                        invϕ,
                        a,
                        b,
                        pb.generators.binf,
                        pb.generators.bsup,
                        method)
                        @show @inbounds A[I[1], I[2]] += swsp(intdata)
        end
        @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end



Nmesh = 300
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, Nmesh; T = T)
basis = P1IntLegendreGenerator(m, T; ordermax = 10)
l = length(basis)
M₀ = zeros(T,l,l)
@time fill_mass_matrix2!(basis, M₀; weight = InvX())