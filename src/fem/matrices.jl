#####################################################################
#                          OVERLAP MATRIX
#####################################################################

function mass_matrix(pb::PolynomialBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_mass_matrix!(pb, A; weight = weight, method = method)
    A
end

function sparse_mass_matrix(pb::PolynomialBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = spzeros(T, size, size)
    fill_mass_matrix!(pb, A; weight = weight, method = method)
    A
end

function mass_matrix(pb::PolynomialBasis, n::Int; method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        mass_matrix(pb; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        mass_matrix(pb; weight = InvX2(), method = default_method(method, InvX2()))
    end
end

function sparse_mass_matrix(pb::PolynomialBasis, n::Int; method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        sparse_mass_matrix(pb; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        sparse_mass_matrix(pb; weight = InvX2(), method = default_method(method, InvX2()))
    end
end

function fill_mass_matrix!(pb::PolynomialBasis, n::Int, A::AbstractArray{<:Real}; method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        fill_mass_matrix!(pb, A; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        fill_mass_matrix!(pb, A; weight = InvX2(), method = default_method(method, InvX2()))
    end
    nothing
end

function fill_mass_matrix!( pb::PolynomialBasis, 
                            A::AbstractMatrix{<:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices, cells_to_generators = pb
    @unpack K = cache
    fill!(A, 0)
    if typeof(weight) == NoWeight
        eldata = getelement(pb, firstindex(mesh),:M)
        _fill_local_matrix!(K, method, weight, eldata, cache.prodMG, pb)
        @inbounds for k ∈ iterators(mesh)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib]
            @views vK = K[Ig, Ig]
            @. vA += vK * invshifts[k][1] / invshifts[1][1]
        end
    else
        @inbounds for k ∈ iterators(mesh)
            eldata = getelement(pb, k,:M)
            _fill_local_matrix!(K, method, weight, eldata, cache.prodMG, pb)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib]
            @views vK = K[Ig, Ig]
            @. vA += vK
        end

    end
    nothing
end


#####################################################################
#                          STIFFNESS MATRIX
#####################################################################

function stiffness_matrix(pb::PolynomialBasis; method::IntegrationMethod = ExactIntegration())
    @unpack size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_stiffness_matrix!(pb, A; method = method)
    A
end

function sparse_stiffness_matrix(pb::PolynomialBasis; method::IntegrationMethod = ExactIntegration())
    @unpack size = pb
    T = eltype(pb)
    A = spzeros(T, size, size)
    fill_stiffness_matrix!(pb, A; method = method)
    A
end

function fill_stiffness_matrix!( pb::PolynomialBasis, 
                                A::AbstractMatrix{<:Real};
                                method::IntegrationMethod = ExactIntegration())
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices, cells_to_generators = pb
    @unpack K = cache
    fill!(A, 0)
    eldata = getelement(pb, firstindex(mesh),:Md)
    _fill_local_matrix!(K, method, NoWeight(), eldata, cache.prodMdG, pb)
    @inbounds for k ∈ iterators(mesh)
        Ib = cells_to_indices[k]
        Ig = cells_to_generators[k]
        @views vA = A[Ib, Ib]
        @views vK = K[Ig, Ig]
        @. vA += vK * shifts[k][1]^2 * invshifts[k][1] / invshifts[1][1]
    end
    nothing
end



#####################################################################
#                          MASS TENSOR
#####################################################################

function mass_tensor(pb::PolynomialBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = default_method(weight))
    T = eltype(pb)
    A = zeros(T, pb.size, pb.size, pb.size)
    fill_mass_tensor!(pb, A; weight = weight, method = method)
    A
end

function mass_tensor(pb::PolynomialBasis, n::Int; method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        mass_tensor(pb; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        mass_tensor(pb; weight = InvX2(), method = default_method(method, InvX2()))
    end
end

function fill_mass_tensor!( pb::PolynomialBasis, 
                            A::AbstractArray{<:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = default_method(weight))
    @unpack mesh = pb
    J = fill(0,3,pb.max_length_intersection[2])
    for I ∈ pb.tensor_fill_indices
        count = find_intersection_indices!(J, pb.indices_cells[I[1]],pb.indices_cells[I[2]], pb.indices_cells[I[3]])
        for c ∈ 1:c
            i = J[1,c]
            j = J[2,c]
            k = J[3,c]
            P = getgenerator(pb, I[1], i)
            Q = getgenerator(pb, I[2], j)
            L = getgenerator(pb, I[3], k)
            ϕ = getshift(pb, I[1], i)
            invϕ = getinvshift(pb, I[1], i)
            (a,b) = getmesh(pb, I[1], i)
            iP = pb.indices_generators[I[1]][i]
            iQ = pb.indices_generators[I[2]][j]
            iL = pb.indices_generators[I[3]][k]
            prod = getprod(pb, iP, iQ, iL)
            intdata = IntegrationData(weight,
                                      (P,Q,L),
                                      prod,
                                      ϕ,
                                      invϕ,
                                      a,
                                      b,
                                      pb.generators.binf,
                                      pb.generators.bsup,
                                      method)
            @inbounds A[I[1], I[2], I[3]] += swsp(intdata)
        end
        @inbounds A[I[3], I[1], I[2]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[2], I[3], I[1]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[2], I[1], I[3]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[3], I[2], I[1]]  = A[I[1], I[2], I[3]]
        @inbounds A[I[1], I[3], I[2]]  = A[I[1], I[2], I[3]]
    end
    nothing
end

function fill_mass_tensor!(pb::PolynomialBasis, 
                           n::Int, 
                           A::Union{AbstractArray{<:Real},Dict{Tuple{Int64, Int64, Int64}, <:Real}};
                           method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        fill_mass_tensor!(pb, A; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        fill_mass_tensor!(pb, A; weight = InvX2(), method = default_method(method, InvX2()))
    end
    nothing
end

function fill_mass_tensor!( pb::PolynomialBasis,
                            A::Dict{Tuple{Int64, Int64, Int64}, <:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices, cells_to_generators = pb
    @unpack T = cache
    fill!(A, 0)
    if typeof(weight) == NoWeight
        eldata = getelement(pb, firstindex(mesh),:T)
        _fill_local_matrix!(T, method, weight, eldata, cache.prodMG, pb)
        @inbounds for k ∈ iterators(mesh)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib]
            @views vK = K[Ig, Ig]
            @. vA += vK
        end
    else
        @inbounds for k ∈ iterators(mesh)
            eldata = getelement(pb, k,:M)
            _fill_local_matrix!(K, method, weight, eldata, cache.prodMG, pb)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib]
            @views vK = K[Ig, Ig]
            @. vA += vK
        end

    end
    nothing
end