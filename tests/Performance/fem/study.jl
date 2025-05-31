using KohnShamResolution
using KohnShamResolution:getgenerator, getshift, getinvshift, IntegrationData, IntegrationMethod, find_intersection_indices, swsp,
        AbstractWeight, NoWeight, InvX
using UnPack

function fill_mass_matrix2!( pb::PolynomialBasis, 
    A::AbstractMatrix{<:Real};
    weight::AbstractWeight = NoWeight(),
    method::IntegrationMethod = ExactIntegration())
    @unpack generators, mesh = pb
    for I ∈ pb.matrix_fill_indices[1:4]
        for (i,j) ∈ find_intersection_indices(pb.indices_cells[I[1],:], pb.indices_cells[I[2],:])
            @show @allocated P = getgenerator(pb, I[1], i)
            @show @allocated Q = getgenerator(pb, I[2], j)
            @show @allocated ϕ = getshift(pb, I[1], i)
            @show @allocated invϕ = getinvshift(pb, I[1], i)
            @show @allocated @views a = mesh[pb.indices_cells[I[1],i]]
            @show @allocated @views b = mesh[pb.indices_cells[I[1],i]+1]
            @show @allocated  intdata = IntegrationData(weight,
                        (P,Q),
                        ϕ,
                        invϕ,
                        a,
                        b,
                        pb.generators.binf,
                        pb.generators.bsup,
                        method)
            @show @allocated  @inbounds A[I[1], I[2]] += swsp(intdata)
        end
        @show @allocated  @inbounds A[I[2],I[1]]  = A[I[1],I[2]]
    end
    nothing
end



Nmesh = 50
T = Float64
Rmin = 0.0
Rmax = 50.0
m = linmesh(Rmin, Rmax, Nmesh; T = T)
basis = P1IntLegendreGenerator(m, T; ordermax = 2)
l = length(basis)

M₀ = zeros(T,l,l)
@code_warntype fill_mass_matrix2!(basis, M₀;weight=InvX())
