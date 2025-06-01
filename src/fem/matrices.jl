#####################################################################
#                          OVERLAP MATRIX
#####################################################################

function mass_matrix(pb::FEMBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_mass_matrix!(pb, A; weight = weight, method = method)
    A
end

function sparse_mass_matrix(pb::FEMBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, size = pb
    T = eltype(pb)
    A = spzeros(T, size, size)
    fill_mass_matrix!(pb, A; weight = weight, method = method)
    A
end

function mass_matrix(pb::FEMBasis, n::Int; method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        mass_matrix(pb; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        mass_matrix(pb; weight = InvX2(), method = default_method(method, InvX2()))
    end
end

function sparse_mass_matrix(pb::FEMBasis, n::Int; method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        sparse_mass_matrix(pb; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        sparse_mass_matrix(pb; weight = InvX2(), method = default_method(method, InvX2()))
    end
end

function fill_mass_matrix!(pb::FEMBasis, n::Int, A::AbstractArray{<:Real}; method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        fill_mass_matrix!(pb, A; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        fill_mass_matrix!(pb, A; weight = InvX2(), method = default_method(method, InvX2()))
    end
    nothing
end

function fill_mass_matrix!( pb::FEMBasis, 
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

function stiffness_matrix(pb::FEMBasis; method::IntegrationMethod = ExactIntegration())
    @unpack size = pb
    T = eltype(pb)
    A = zeros(T, size, size)
    fill_stiffness_matrix!(pb, A; method = method)
    A
end

function sparse_stiffness_matrix(pb::FEMBasis; method::IntegrationMethod = ExactIntegration())
    @unpack size = pb
    T = eltype(pb)
    A = spzeros(T, size, size)
    fill_stiffness_matrix!(pb, A; method = method)
    A
end

function fill_stiffness_matrix!( pb::FEMBasis, 
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

function mass_tensor(pb::FEMBasis; weight::AbstractWeight = NoWeight(), method::IntegrationMethod = default_method(weight))
    T = eltype(pb)
    A = zeros(T, pb.size, pb.size, pb.size)
    fill_mass_tensor!(pb, A; weight = weight, method = method)
    A
end

function mass_tensor(pb::FEMBasis, n::Int; method::IntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        mass_tensor(pb; weight = InvX(), method = default_method(method, InvX()))
    elseif n == -2
        mass_tensor(pb; weight = InvX2(), method = default_method(method, InvX2()))
    end
end

function fill_mass_tensor!( pb::FEMBasis, 
                            A::AbstractArray{<:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices, cells_to_generators = pb
    @unpack T = cache
    fill!(A, 0)
    if typeof(weight) == NoWeight
        eldata = getelement(pb, firstindex(mesh),:T)
        _fill_local_matrix!(T, method, weight, eldata, cache.prodTG, pb)
        @inbounds for k ∈ iterators(mesh)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib, Ib]
            @views vT = T[Ig, Ig, Ig]
            @. vA += vT * invshifts[k][1] / invshifts[1][1]
        end
    else
        @inbounds for k ∈ iterators(mesh)
            eldata = getelement(pb, k,:T)
            _fill_local_matrix!(T, method, weight, eldata, cache.prodTG, pb)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib, Ib]
            @views vT = T[Ig, Ig, Ig]
            @. vA += vT
        end
    end
    nothing
end

function fill_mass_tensor!(pb::FEMBasis, 
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

function fill_mass_tensor!( pb::FEMBasis,
                            A::Dict{Tuple{Int64, Int64, Int64}, <:Real};
                            weight::AbstractWeight = NoWeight(),
                            method::IntegrationMethod = default_method(weight))
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices, cells_to_generators = pb
    @unpack T = cache
    empty!(A)  # Réinitialise proprement le dictionnaire

    if typeof(weight) == NoWeight
        eldata = getelement(pb, firstindex(mesh), :T)
        _fill_local_matrix!(T, method, weight, eldata, cache.prodTG, pb)
        @inbounds for k ∈ iterators(mesh)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            scaling = invshifts[k][1] / invshifts[1][1]
            for (ii, i) in enumerate(Ib), (jj, j) in enumerate(Ib), (kk, k3) in enumerate(Ib)
                key = (i, j, k3)
                val = T[Ig[ii], Ig[jj], Ig[kk]] * scaling
                if haskey(A, key)
                    A[key] += val
                else
                    A[key] = val
                end
            end
        end
    else
        @inbounds for k ∈ iterators(mesh)
            eldata = getelement(pb, k, :T)
            _fill_local_matrix!(T, method, weight, eldata, cache.prodTG, pb)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            for (ii, i) in enumerate(Ib), (jj, j) in enumerate(Ib), (kk, k3) in enumerate(Ib)
                key = (i, j, k3)
                val = T[Ig[ii], Ig[jj], Ig[kk]]
                if haskey(A, key)
                    A[key] += val
                else
                    A[key] = val
                end
            end
        end
    end
    return nothing
end