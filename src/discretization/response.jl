# ===================================================================
#                 LINEAR RESPONSE (LDA / LSDA ONLY)
# ===================================================================
# Density matrices and Hamiltonian variations are (Nₕ, Nₕ, nspin) arrays here,
# also for nspin = 1, where the SCF density matrix is an (Nₕ, Nₕ) matrix.
_spin(A::AbstractMatrix, σ::Int) = A
_spin(A::AbstractArray{<:Real, 3}, σ::Int) = view(A, :, :, σ)

"""
    ResponseWorkspace(discretization)

Buffers at the quadrature points of one cell for the linear response: the
density `ρ`, its variation `δρ`, the exchange-correlation kernel and a cache for
its correlation part (Libxc `v2rho2` layout), the Hartree response matrix, and
the kernel at the quadrature points of every cell (`kernels[k]`), filled by
[`prepare_kernel!`](@ref).
"""
struct ResponseWorkspace{T, Tρ, Tk}
    ρ::Tρ
    δρ::Tρ
    kernel::Tk
    kcache::Tk
    VH::Matrix{T}
    kernels::Vector{Tk}
    function ResponseWorkspace(discretization::KSEDiscretization{T}) where T
        @unpack nspin, cache = discretization
        neval = size(cache.excw.ρ_buf, ndims(cache.excw.ρ_buf))
        ρ = flexible_zeros(T, nspin, (neval,))
        δρ = flexible_zeros(T, nspin, (neval,))
        nk = nspin == 1 ? 1 : 3
        kernel = flexible_zeros(T, nk, (neval,))
        kcache = flexible_zeros(T, nk, (neval,))
        VH = zeros(T, discretization.Nₕ, discretization.Nₕ)
        kernels = [copy(kernel) for _ in cellrange(discretization.basis.mesh)]
        new{T, typeof(ρ), typeof(kernel)}(ρ, δρ, kernel, kcache, VH, kernels)
    end
end

"""
    prepare_kernel!(ws, discretization, model, D)

Exchange-correlation kernel `f_xc(ρ)` of the density matrix `D` at the quadrature
points of every cell, stored in `ws.kernels`: it does not depend on the density
variation, so [`assemble_exc_response!`](@ref) can reuse it for any `δD`.
"""
function prepare_kernel!(ws::ResponseWorkspace, discretization::KSEDiscretization,
                         model::KSEModel, D::AbstractArray{<:Real})
    @unpack basis, nspin, fem_integration_method = discretization
    @unpack x = fem_integration_method
    @unpack ρ, kcache = ws
    has_gga(model) && error("prepare_kernel!: GGA is not implemented.")
    X = similar(x)
    for k in cellrange(basis.mesh)
        eldata = getelement(basis, k, :M)
        @. X = eldata.invϕ[1] * x + eldata.invϕ[2]
        if nspin == 1
            optimized_eval_density!(ρ, discretization, _spin(D, 1), X)
        else
            for σ in 1:2
                @views optimized_eval_density!(ρ[σ, :], discretization, D[:, :, σ], X)
            end
        end
        evaluate_v2rho2!(model; rho = ρ, v2rho2 = ws.kernels[k], cache = kcache)
    end
    ws
end

"""
    assemble_exc_response!(δV, ws, discretization, model, δD)

Exchange-correlation response matrices `δV[:, :, σ] = ∫ χᵢ χⱼ Σσ' f_xc^{σσ'}(ρ) δρσ'`
for a density-matrix variation `δD`, with the kernel precomputed by
[`prepare_kernel!`](@ref), assembled with the same quadrature as
[`assemble_exc!`](@ref). LDA and LSDA only.
"""
function assemble_exc_response!(δV::AbstractArray{<:Real}, ws::ResponseWorkspace,
                                discretization::KSEDiscretization, model::KSEModel,
                                δD::AbstractArray{<:Real})
    @unpack basis, nspin, fem_integration_method = discretization
    @unpack δρ, kernels = ws
    mesh = basis.mesh
    if nspin == 1
        δD1 = _spin(δD, 1)
        function _weight!(Y::AbstractVector, X::AbstractVector)
            kern = kernels[findindex(mesh, first(X))]
            optimized_eval_density!(δρ, discretization, δD1, X)
            @. Y = kern * δρ
        end
        weight = FunWeight(_weight!; is_inplace = true, is_vectorized = true)
        @views fill_mass_matrix!(basis, δV[:, :, 1]; weight = weight,
                                 method = fem_integration_method)
    else
        function _weight_up!(Y::AbstractVector, X::AbstractVector)
            kern = kernels[findindex(mesh, first(X))]
            for σ in 1:2
                @views optimized_eval_density!(δρ[σ, :], discretization, δD[:, :, σ], X)
            end
            @views @. Y = kern[1, :] * δρ[1, :] + kern[2, :] * δρ[2, :]
        end
        function _weight_down!(Y::AbstractVector, X::AbstractVector)
            kern = kernels[findindex(mesh, first(X))]
            for σ in 1:2
                @views optimized_eval_density!(δρ[σ, :], discretization, δD[:, :, σ], X)
            end
            @views @. Y = kern[2, :] * δρ[1, :] + kern[3, :] * δρ[2, :]
        end
        wup = FunWeight(_weight_up!; is_inplace = true, is_vectorized = true)
        wdown = FunWeight(_weight_down!; is_inplace = true, is_vectorized = true)
        @views fill_mass_matrix!(basis, δV[:, :, 1]; weight = wup, method = fem_integration_method)
        @views fill_mass_matrix!(basis, δV[:, :, 2]; weight = wdown, method = fem_integration_method)
    end
    for σ in 1:nspin
        @views δV[:, :, σ] .= (δV[:, :, σ] .+ δV[:, :, σ]') ./ 2
    end
    δV
end

"""
    apply_kernel!(δH, ws, discretization, model, D, δD; prepare = true)

Variation `δH[:, :, σ]` of the mean-field Hamiltonian for a density-matrix
variation `δD` around `D`: Hartree response (with `N` fixed, as in
[`assemble_hamiltonian!`](@ref)) plus exchange-correlation response. With
`prepare = false`, the kernel of `D` must already be in `ws` (see
[`prepare_kernel!`](@ref)); this is how repeated applications at the same `D`
avoid recomputing it.
"""
function apply_kernel!(δH::AbstractArray{<:Real}, ws::ResponseWorkspace,
                       discretization::KSEDiscretization, model::KSEModel,
                       D::AbstractArray{<:Real}, δD::AbstractArray{<:Real};
                       prepare::Bool = true)
    @unpack nspin = discretization
    T = eltype(δH)
    fill!(δH, zero(T))
    if has_exchcorr(model)
        prepare && prepare_kernel!(ws, discretization, model, D)
        assemble_exc_response!(δH, ws, discretization, model, δD)
    end
    if !iszero(model.hartree)
        δVH = ws.VH
        assemble_hartree!(δVH, discretization, nspin == 1 ? _spin(δD, 1) : δD, zero(T),
                          model.hartree)
        for σ in 1:nspin
            @views δH[:, :, σ] .+= δVH
        end
    end
    δH
end

"""
    ChannelData

Complete eigendecomposition of one channel `(l, σ)` of the Hamiltonian: all the
orbitals `U` (columns, `Uᵀ M₀ U = I`) and energies `ϵ`, the indices `occ` of the
occupied orbitals with their occupations `nocc` (including the `2l+1` weight),
and the sum-over-states coefficients `C[a, j] = nocc[j]/(ϵ[occ[j]] - ϵ[a])`
(zero when `a` is occupied with the same occupation, since such pairs cancel).
"""
mutable struct ChannelData{T <: Real}
    l::Int
    σ::Int
    U::Matrix{T}
    ϵ::Vector{T}
    occ::Vector{Int}
    nocc::Vector{T}
    C::Matrix{T}
end

"""
    response_coefficients!(ch::ChannelData)

Fill the sum-over-states coefficients `ch.C` from `ch.ϵ`, `ch.occ`, `ch.nocc`.
"""
function response_coefficients!(ch::ChannelData{T}) where T
    nb = length(ch.ϵ)
    ch.C = zeros(T, nb, length(ch.occ))
    nall = zeros(T, nb)
    nall[ch.occ] .= ch.nocc
    for (j, i) in enumerate(ch.occ), a in 1:nb
        (a == i || nall[a] == nall[i]) && continue
        ch.C[a, j] = ch.nocc[j] / (ch.ϵ[i] - ch.ϵ[a])
    end
    ch
end

"""
    apply_chi0!(δD, channels, δH)

Independent-particle response: the density-matrix variation `δD[:, :, σ]` of the
occupied orbitals of all `channels` under a Hamiltonian variation `δH[:, :, σ]`,
at fixed occupations, `δD = Σ_j Σ_a C[a,j] ⟨u_a|δH|u_j⟩ (u_a u_jᵀ + u_j u_aᵀ)`.
"""
function apply_chi0!(δD::AbstractArray{<:Real}, channels::Vector{<:ChannelData},
                     δH::AbstractArray{<:Real})
    fill!(δD, zero(eltype(δD)))
    for ch in channels
        isempty(ch.occ) && continue
        @views Uocc = ch.U[:, ch.occ]
        @views G = ch.U' * (δH[:, :, ch.σ] * Uocc)
        G .*= ch.C
        Y = ch.U * G
        YU = Y * Uocc'
        @views δD[:, :, ch.σ] .+= YU .+ YU'
    end
    δD
end

"""
    diagonalize_channels!(channels, discretization, H, δH = nothing)

Full eigendecomposition of every channel `(l, σ)` of `H` (optionally `H + δH`,
with `δH[:, :, σ]` independent of `l`), stored in `channels` (orbitals, energies;
occupations and coefficients are left untouched).
"""
function diagonalize_channels!(channels::Vector{<:ChannelData},
                               discretization::KSEDiscretization,
                               H::AbstractArray{<:Real}, δH = nothing)
    @unpack S = discretization.femops
    @threads for ch in channels
        @views Hc = H[:, :, ch.l + 1, ch.σ]
        A = δH === nothing ? S * Hc * S : S * (Hc .+ _spin(δH, ch.σ)) * S
        λ, V = eigen(Symmetric(A))
        ch.ϵ .= λ
        mul!(ch.U, S, V)
    end
    channels
end

"""
    density_matrix!(D, channels)

Density matrix `D[:, :, σ] = Σ_j nocc[j] u_j u_jᵀ` of the occupied orbitals of
`channels`.
"""
function density_matrix!(D::AbstractArray{<:Real}, channels::Vector{<:ChannelData})
    fill!(D, zero(eltype(D)))
    for ch in channels, (j, i) in enumerate(ch.occ)
        @views u = ch.U[:, i]
        @views _spin(D, ch.σ) .+= ch.nocc[j] .* (u * u')
    end
    D
end

"""
    empty_channels(discretization) -> Vector{ChannelData}

One `ChannelData` per channel `(l, σ)`, with no occupied orbital.
"""
function empty_channels(discretization::KSEDiscretization{T}) where T
    @unpack lₕ, nspin, Nₕ = discretization
    [ChannelData{T}(l, σ, zeros(T, Nₕ, Nₕ), zeros(T, Nₕ), Int[], T[], zeros(T, Nₕ, 0))
     for σ in 1:nspin for l in 0:lₕ]
end
