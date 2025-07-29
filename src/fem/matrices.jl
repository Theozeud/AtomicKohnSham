#--------------------------------------------------------------------
#                           ELEMENT DATA STRUCTURE
#--------------------------------------------------------------------
struct ElementData{T <: Real}
    ϕ::Tuple{T, T}
    invϕ::Tuple{T, T}
    a::T
    b::T
    binf::T
    bsup::T
    s::Symbol

    function ElementData(ϕ::Tuple{<:Real, <:Real},
            invϕ::Tuple{<:Real, <:Real},
            a::Real,
            b::Real,
            binf::Real,
            bsup::Real,
            s::Symbol)
        new{typeof(a)}(ϕ,
            invϕ,
            a,
            b,
            binf,
            bsup,
            s)
    end
end

function ElementData(eldata::ElementData;
        ϕ::Tuple{<:Real, <:Real} = eldata.ϕ,
        invϕ::Tuple{<:Real, <:Real} = eldata.invϕ,
        a::Real = eldata.a,
        b::Real = eldata.b,
        binf::Real = eldata.binf,
        bsup::Real = eldata.bsup,
        s::Symbol = eldata.s)
    ElementData(ϕ,
        invϕ,
        a,
        b,
        binf,
        bsup,
        s)
end

function getelement(basis::FEMBasis, k::Int, s::Symbol)
    @unpack generators, mesh, shifts, invshifts = basis
    ElementData(shifts[k],
        invshifts[k],
        mesh[k],
        mesh[k + 1],
        generators.binf,
        generators.bsup,
        s)
end

#--------------------------------------------------------------------
#                          OVERLAP MATRIX
#--------------------------------------------------------------------

function mass_matrix(basis::FEMBasis;
        weight::AbstractWeight = NoWeight(),
        method::FEMIntegrationMethod = default_method(basis, weight))
    @unpack generators, mesh, size = basis
    T = eltype(basis)
    A = zeros(T, size, size)
    fill_mass_matrix!(basis, A; weight = weight, method = method)
    A
end

function sparse_mass_matrix(basis::FEMBasis;
        weight::AbstractWeight = NoWeight(),
        method::FEMIntegrationMethod = default_method(basis, weight))
    @unpack generators, mesh, size = basis
    T = eltype(basis)
    A = spzeros(T, size, size)
    fill_mass_matrix!(basis, A; weight = weight, method = method)
    A
end

function mass_matrix(basis::FEMBasis,
        n::Int;
        method::FEMIntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        mass_matrix(basis; weight = InvX(), method = default_method(basis, method, InvX()))
    elseif n == -2
        mass_matrix(
            basis; weight = InvX2(), method = default_method(basis, method, InvX2()))
    end
end

function sparse_mass_matrix(basis::FEMBasis,
        n::Int;
        method::FEMIntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        sparse_mass_matrix(
            basis; weight = InvX(), method = default_method(basis, method, InvX()))
    elseif n == -2
        sparse_mass_matrix(
            basis; weight = InvX2(), method = default_method(basis, method, InvX2()))
    end
end

function fill_mass_matrix!(basis::FEMBasis, n::Int, A::AbstractArray{<:Real};
        method::FEMIntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        fill_mass_matrix!(
            basis, A; weight = InvX(), method = default_method(basis, method, InvX()))
    elseif n == -2
        fill_mass_matrix!(
            basis, A; weight = InvX2(), method = default_method(basis, method, InvX2()))
    end
    nothing
end

function fill_mass_matrix!(basis::FEMBasis,
        A::AbstractMatrix{<:Real};
        weight::AbstractWeight = NoWeight(),
        method::FEMIntegrationMethod = default_method(basis, weight))
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices,
    cells_to_generators = basis
    @unpack K = cache
    fill!(A, 0)
    if typeof(weight) == NoWeight
        eldata = getelement(basis, firstindex(mesh), :M)
        _fill_local_matrix!(K, method, weight, eldata, cache.prodMG, basis)
        @inbounds for k in cellrange(mesh)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib]
            @views vK = K[Ig, Ig]
            @. vA += vK * invshifts[k][1] / invshifts[1][1]
        end
    else
        @inbounds for k in cellrange(mesh)
            eldata = getelement(basis, k, :M)
            _fill_local_matrix!(K, method, weight, eldata, cache.prodMG, basis)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib]
            @views vK = K[Ig, Ig]
            @. vA += vK
        end
    end
    nothing
end

#--------------------------------------------------------------------
#                          STIFFNESS MATRIX
#--------------------------------------------------------------------
function stiffness_matrix(
        basis::FEMBasis; method::FEMIntegrationMethod = ExactIntegration())
    @unpack size = basis
    T = eltype(basis)
    A = zeros(T, size, size)
    fill_stiffness_matrix!(basis, A; method = method)
    A
end

function sparse_stiffness_matrix(
        basis::FEMBasis; method::FEMIntegrationMethod = ExactIntegration())
    @unpack size = basis
    T = eltype(basis)
    A = spzeros(T, size, size)
    fill_stiffness_matrix!(basis, A; method = method)
    A
end

function fill_stiffness_matrix!(basis::FEMBasis,
        A::AbstractMatrix{<:Real};
        method::FEMIntegrationMethod = ExactIntegration())
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices,
    cells_to_generators = basis
    @unpack K = cache
    fill!(A, 0)
    eldata = getelement(basis, firstindex(mesh), :Md)
    _fill_local_matrix!(K, method, NoWeight(), eldata, cache.prodMdG, basis)
    @inbounds for k in cellrange(mesh)
        Ib = cells_to_indices[k]
        Ig = cells_to_generators[k]
        @views vA = A[Ib, Ib]
        @views vK = K[Ig, Ig]
        @. vA += vK * shifts[k][1]^2 * invshifts[k][1] / invshifts[1][1]
    end
    nothing
end

#--------------------------------------------------------------------
#                          MASS TENSOR
#--------------------------------------------------------------------

function mass_tensor(basis::FEMBasis; weight::AbstractWeight = NoWeight(),
        method::FEMIntegrationMethod = default_method(basis, weight))
    T = eltype(basis)
    A = zeros(T, basis.size, basis.size, basis.size)
    fill_mass_tensor!(basis, A; weight = weight, method = method)
    A
end

function mass_tensor(
        basis::FEMBasis, n::Int; method::FEMIntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        mass_tensor(basis; weight = InvX(), method = default_method(basis, method, InvX()))
    elseif n == -2
        mass_tensor(
            basis; weight = InvX2(), method = default_method(basis, method, InvX2()))
    end
end

function fill_mass_tensor!(basis::FEMBasis,
        A::AbstractArray{<:Real};
        weight::AbstractWeight = NoWeight(),
        method::FEMIntegrationMethod = default_method(basis, weight))
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices,
    cells_to_generators = basis
    @unpack T = cache
    fill!(A, 0)
    if typeof(weight) == NoWeight
        eldata = getelement(basis, firstindex(mesh), :T)
        _fill_local_matrix!(T, method, weight, eldata, cache.prodTG, basis)
        @inbounds for k in cellrange(mesh)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib, Ib]
            @views vT = T[Ig, Ig, Ig]
            @. vA += vT * invshifts[k][1] / invshifts[1][1]
        end
    else
        @inbounds for k in cellrange(mesh)
            eldata = getelement(basis, k, :T)
            _fill_local_matrix!(T, method, weight, eldata, cache.prodTG, basis)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            @views vA = A[Ib, Ib, Ib]
            @views vT = T[Ig, Ig, Ig]
            @. vA += vT
        end
    end
    nothing
end

function fill_mass_tensor!(basis::FEMBasis,
        n::Int,
        A::Union{AbstractArray{<:Real}, Dict{Tuple{Int64, Int64, Int64}, <:Real}};
        method::FEMIntegrationMethod = NoSelectedMethod())
    @assert n == -1 || n == -2 "n must be -1 or -2"
    if n == -1
        fill_mass_tensor!(
            basis, A; weight = InvX(), method = default_method(basis, method, InvX()))
    elseif n == -2
        fill_mass_tensor!(
            basis, A; weight = InvX2(), method = default_method(basis, method, InvX2()))
    end
    nothing
end

function fill_mass_tensor!(basis::FEMBasis,
        A::Dict{Tuple{Int64, Int64, Int64}, <:Real};
        weight::AbstractWeight = NoWeight(),
        method::FEMIntegrationMethod = default_method(basis, weight))
    @unpack generators, mesh, cache, shifts, invshifts, cells_to_indices,
    cells_to_generators = basis
    @unpack T = cache
    empty!(A)  # Réinitialise proprement le dictionnaire

    if typeof(weight) == NoWeight
        eldata = getelement(basis, firstindex(mesh), :T)
        _fill_local_matrix!(T, method, weight, eldata, cache.prodTG, basis)
        @inbounds for k in cellrange(mesh)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            scaling = invshifts[k][1] / invshifts[1][1]
            for (ii, i) in enumerate(Ib), (jj, j) in enumerate(Ib),
                (kk, k3) in enumerate(Ib)
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
        @inbounds for k in cellrange(mesh)
            eldata = getelement(basis, k, :T)
            _fill_local_matrix!(T, method, weight, eldata, cache.prodTG, basis)
            Ib = cells_to_indices[k]
            Ig = cells_to_generators[k]
            for (ii, i) in enumerate(Ib), (jj, j) in enumerate(Ib),
                (kk, k3) in enumerate(Ib)
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
