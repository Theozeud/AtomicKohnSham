# ===================================================================
#                          Kinetic Matrix
# ===================================================================
"""
Assemble the kinetic operators for all angular momentum channels.
"""
function assemble_kinetic!(discretization::KSEDiscretization)
    @unpack A, M₋₂ = discretization.femops
    @unpack Kin = discretization.ksham
    for l in 0:discretization.lₕ
        @. Kin[l + 1] = (A + l*(l+1)*M₋₂) /2
    end
    nothing
end

# ===================================================================
#                          Coulomb Matrix
# ===================================================================
"""
Assemble the nuclear Coulomb potential operator for charge `Z`.
"""
function assemble_coulomb!(discretization::KSEDiscretization, Z::Real)
    @unpack M₋₁ = discretization.femops
    @unpack Coulomb = discretization.ksham
    Coulomb .= - Z .* M₋₁
    nothing
end

# ===================================================================
#                          Hartree Matrix
# ===================================================================
"""
Assemble the Hartree potential function from the current density `D`.
"""
function assemble_hartree_pot!(discretization::KSEDiscretization,
                              D::AbstractArray{<:Real}; coeff::Real = true)
    @unpack Rmax, cache, nspin = discretization
    @unpack A, F = discretization.femops
    @unpack Hartree = discretization.ksham
    @unpack B, W = cache.hartw
    # HartreeWorkspace allocates B/W with length 0 when model.hartree == 0
    # (see discretization.jl); A\B below would otherwise throw a
    # DimensionMismatch. W is already empty in that case, so there's nothing
    # left to do.
    isempty(B) && return nothing
    if nspin == 1
        tensor_matrix_dict!(B, D, F)
        W .= A\B
    else
        @views DUP = D[:, :, 1]
        @views DDOWN = D[:, :, 2]
        tensor_matrix_dict!(B, DUP, DDOWN, F)
        W .= A\B
    end
    W .*= coeff
    nothing
end


"""
Assemble the Hartree potential operator from the current density `D`.
"""
function assemble_hartree!(discretization::KSEDiscretization,
                          D::AbstractArray{<:Real}, N::Real, coeff::Real = true)
    @unpack Rmax, cache, nspin = discretization
    @unpack A, F, M₀ = discretization.femops
    @unpack Hartree = discretization.ksham
    @unpack FW, B, W = cache.hartw
    if nspin == 1
        tensor_matrix_dict!(B, D, F)
        W .= A\B
    else
        @views DUP = D[:, :, 1]
        @views DDOWN = D[:, :, 2]
        tensor_matrix_dict!(B, DUP, DDOWN, F)
        W .= A\B
    end
    tensor_vector_dict!(FW, W, F)
    @. Hartree = FW + N/Rmax * M₀
    @. Hartree .*= coeff
    Hartree .= (Hartree .+ Hartree') ./ 2
    nothing
end

# ===================================================================
#                   Exchange Correlation Matrix
# ===================================================================
"""
Evaluate the radial density at quadrature points from a local FEM stencil.
"""
function optimized_eval_density!(ρ::AbstractVector{<:Real},
                                discretization::KSEDiscretization,
                                D::AbstractMatrix{<:Real},
                                X::AbstractVector{<:Real})
    @unpack basis, fem_integration_method= discretization
    Qgenx = fem_integration_method.Qgenx
    fill!(ρ, 0)
    idxmesh = findindex(basis.mesh, first(X))
    Ib = basis.cells_to_indices[idxmesh]
    Ig = basis.cells_to_generators[idxmesh]
    a = Int(sqrt(size(Qgenx, 1)))
    Qgenxreshape = reshape(Qgenx, a, a, size(Qgenx, 2))
    # Ib/Ig are Vector{Int} permutations of a contiguous block (not ranges,
    # see cells_to_indices/cells_to_generators), so D[Ib, Ib] and
    # Qgenxreshape[Ig, Ig, :] are non-strided views: @tensor can't contract
    # them via BLAS and would silently materialize dense copies every call
    # (~neval*a^2 elements) instead. A scalar loop indexes them directly, with
    # no temporary allocation.
    n = length(Ib)
    neval = length(ρ)
    @inbounds for j in 1:n
        Ibj = Ib[j]
        Igj = Ig[j]
        for i in 1:n
            Dij = D[Ib[i], Ibj]
            Igi = Ig[i]
            for k in 1:neval
                ρ[k] += Dij * Qgenxreshape[Igi, Igj, k]
            end
        end
    end
    @.ρ /= X^2
    @. ρ *= 1/4π
    ρ
end

"""
Evaluate the radial density gradient ρ'(x) at quadrature points from a local
FEM stencil (GGA). Requires `ρ` already evaluated at the same points `X` (see
[`optimized_eval_density!`](@ref)), since `ρ'(r) = q'(r)/(4πr²) - 2ρ(r)/r`
needs `ρ` -- avoids a redundant contraction of the same density matrix.
"""
function optimized_eval_density_gradient!(dρ::AbstractVector{<:Real},
                                discretization::KSEDiscretization,
                                D::AbstractMatrix{<:Real},
                                ρ::AbstractVector{<:Real},
                                X::AbstractVector{<:Real})
    @unpack basis, fem_integration_method = discretization
    Qmixedgenx = fem_integration_method.Qmixedgenx
    fill!(dρ, 0)
    idxmesh = findindex(basis.mesh, first(X))
    Ib = basis.cells_to_indices[idxmesh]
    Ig = basis.cells_to_generators[idxmesh]
    a = Int(sqrt(size(Qmixedgenx, 1)))
    Qmixedgenxreshape = reshape(Qmixedgenx, a, a, size(Qmixedgenx, 2))
    # See optimized_eval_density! for why this is a scalar loop, not @tensor.
    n = length(Ib)
    neval = length(dρ)
    @inbounds for j in 1:n
        Ibj = Ib[j]
        Igj = Ig[j]
        for i in 1:n
            Dij = D[Ib[i], Ibj]
            Igi = Ig[i]
            for k in 1:neval
                dρ[k] += Dij * Qmixedgenxreshape[Igi, Igj, k]
            end
        end
    end
    # Qmixedgenx holds φᵢ'(reference)·φⱼ(reference): unlike Qgenx (used in
    # optimized_eval_density!), a single derivative factor needs the
    # physical<-reference chain-rule slope basis.shifts[idxmesh][1] (same ϕ[1]
    # factor as evaluate_deriv! in fem/basis.jl) before it's a physical q'(r).
    ϕ1 = basis.shifts[idxmesh][1]
    @. dρ = 2*ϕ1*dρ/(4π*X^2) - 2ρ/X
    dρ
end

"""
Assemble the exchange–correlation potential operators from the current density `D`.

For a GGA model (`has_gga(model)`), the weak-form derivation in `dev/gga/PLAN.md`
gives, per spin channel σ ∈ {↑,↓}:

    Vxcσ = M[vrhoσ - 2Kσ/r] + Mmix[Kσ] + Mmix[Kσ]ᵀ

with `K↑ = 2vσ↑↑·ρ↑' + vσ↑↓·ρ↓'` (and ↑↔↓ swapped for `K↓`); for `nspin == 1`
there is only the "self" term, `K = 2vσ·ρ'`.
"""
function assemble_exc!(discretization::KSEDiscretization, model::KSEModel,
                      D::AbstractArray{<:Real})
    @unpack basis, nspin, fem_integration_method, cache = discretization
    @unpack ρ_buf, vρ_buf, vρ_buf2, dρ_buf, σ_buf, vσ_buf, vσ_buf2 = cache.excw
    @unpack VxcUP, VxcDOWN, VxcMix = discretization.ksham
    if has_gga(model)
        if nspin == 1
            function _weight_gga!(Y::AbstractVector, X::AbstractVector)
                optimized_eval_density!(ρ_buf, discretization, D, X)
                optimized_eval_density_gradient!(dρ_buf, discretization, D, ρ_buf, X)
                @. σ_buf = dρ_buf^2
                evaluate_vrho_vsigma!(model; rho = ρ_buf, sigma = σ_buf,
                                     vrho = Y, vsigma = vσ_buf,
                                     cache_vrho = vρ_buf, cache_vsigma = vσ_buf2)
                @. Y = Y - 4 * vσ_buf * dρ_buf / X
            end
            function _weightmix_gga!(Y::AbstractVector, X::AbstractVector)
                optimized_eval_density!(ρ_buf, discretization, D, X)
                optimized_eval_density_gradient!(dρ_buf, discretization, D, ρ_buf, X)
                @. σ_buf = dρ_buf^2
                evaluate_vrho_vsigma!(model; rho = ρ_buf, sigma = σ_buf,
                                     vrho = vρ_buf2, vsigma = vσ_buf,
                                     cache_vrho = vρ_buf, cache_vsigma = vσ_buf2)
                @. Y = 2 * vσ_buf * dρ_buf
            end
            weight = FunWeight(_weight_gga!; is_inplace = true, is_vectorized = true)
            weightmix = FunWeight(_weightmix_gga!; is_inplace = true, is_vectorized = true)
            fill!(VxcUP, 0)
            fill_mass_matrix!(basis, VxcUP; weight = weight, method = fem_integration_method)
            fill_mixed_mass_matrix!(
                basis, VxcMix; weight = weightmix, method = fem_integration_method)
            @. VxcUP += VxcMix + VxcMix'
            VxcUP .= (VxcUP .+ VxcUP') ./ 2
        else
            @views DUP = D[:, :, 1]
            @views DDOWN = D[:, :, 2]
            @views ρ_buf_up = ρ_buf[1, :]
            @views ρ_buf_down = ρ_buf[2, :]
            @views dρ_buf_up = dρ_buf[1, :]
            @views dρ_buf_down = dρ_buf[2, :]
            @views σ_uu = σ_buf[1, :]
            @views σ_ud = σ_buf[2, :]
            @views σ_dd = σ_buf[3, :]
            @views vσ_uu = vσ_buf[1, :]
            @views vσ_ud = vσ_buf[2, :]
            @views vσ_dd = vσ_buf[3, :]
            function _eval_common!(X::AbstractVector)
                optimized_eval_density!(ρ_buf_up, discretization, DUP, X)
                optimized_eval_density!(ρ_buf_down, discretization, DDOWN, X)
                optimized_eval_density_gradient!(dρ_buf_up, discretization, DUP, ρ_buf_up, X)
                optimized_eval_density_gradient!(
                    dρ_buf_down, discretization, DDOWN, ρ_buf_down, X)
                @. σ_uu = dρ_buf_up^2
                @. σ_ud = dρ_buf_up * dρ_buf_down
                @. σ_dd = dρ_buf_down^2
                evaluate_vrho_vsigma!(model; rho = ρ_buf, sigma = σ_buf,
                                     vrho = vρ_buf2, vsigma = vσ_buf,
                                     cache_vrho = vρ_buf, cache_vsigma = vσ_buf2)
            end
            function _weight_up_gga!(Y::AbstractVector, X::AbstractVector)
                _eval_common!(X)
                @views vrho_up = vρ_buf2[1, :]
                @. Y = vrho_up - 4 * vσ_uu * dρ_buf_up / X - 2 * vσ_ud * dρ_buf_down / X
            end
            function _weightmix_up_gga!(Y::AbstractVector, X::AbstractVector)
                _eval_common!(X)
                @. Y = 2 * vσ_uu * dρ_buf_up + vσ_ud * dρ_buf_down
            end
            function _weight_down_gga!(Y::AbstractVector, X::AbstractVector)
                _eval_common!(X)
                @views vrho_down = vρ_buf2[2, :]
                @. Y = vrho_down - 4 * vσ_dd * dρ_buf_down / X - 2 * vσ_ud * dρ_buf_up / X
            end
            function _weightmix_down_gga!(Y::AbstractVector, X::AbstractVector)
                _eval_common!(X)
                @. Y = 2 * vσ_dd * dρ_buf_down + vσ_ud * dρ_buf_up
            end
            weightUP = FunWeight(_weight_up_gga!; is_inplace = true, is_vectorized = true)
            weightDOWN = FunWeight(_weight_down_gga!; is_inplace = true, is_vectorized = true)
            weightmixUP = FunWeight(_weightmix_up_gga!; is_inplace = true, is_vectorized = true)
            weightmixDOWN = FunWeight(_weightmix_down_gga!; is_inplace = true, is_vectorized = true)
            fill!(VxcUP, 0)
            fill!(VxcDOWN, 0)
            fill_mass_matrix!(basis, VxcUP; weight = weightUP, method = fem_integration_method)
            fill_mixed_mass_matrix!(
                basis, VxcMix; weight = weightmixUP, method = fem_integration_method)
            @. VxcUP += VxcMix + VxcMix'
            fill_mass_matrix!(
                basis, VxcDOWN; weight = weightDOWN, method = fem_integration_method)
            fill_mixed_mass_matrix!(
                basis, VxcMix; weight = weightmixDOWN, method = fem_integration_method)
            @. VxcDOWN += VxcMix + VxcMix'
            VxcUP .= (VxcUP .+ VxcUP') ./ 2
            VxcDOWN .= (VxcDOWN .+ VxcDOWN') ./ 2
        end
        return nothing
    end
    if nspin == 1
        function _weight!(Y::AbstractVector, X::AbstractVector)
            optimized_eval_density!(ρ_buf, discretization, D, X)
            evaluate_vrho!(model; vrho = Y, rho = ρ_buf, cache = vρ_buf)
        end
        weight = FunWeight(_weight!; is_inplace = true, is_vectorized = true)
        fill!(VxcUP, 0)
        fill_mass_matrix!(basis, VxcUP; weight = weight, method = fem_integration_method)
        VxcUP .= (VxcUP .+ VxcUP') ./ 2
    else
        @views DUP = D[:, :, 1]
        @views DDOWN = D[:, :, 2]
        function _weight_up!(Y::AbstractVector, X::AbstractVector)
            @views ρ_buf_up = ρ_buf[1, :]
            @views ρ_buf_down = ρ_buf[2, :]
            optimized_eval_density!(ρ_buf_up, discretization, DUP, X)
            optimized_eval_density!(ρ_buf_down, discretization, DDOWN, X)
            evaluate_vrho!(model; vrho = vρ_buf2, rho = ρ_buf, cache = vρ_buf)
            @views vρ_buf2up = vρ_buf2[1, :]
            Y .= vρ_buf2up
        end

        function _weight_down!(Y::AbstractVector, X::AbstractVector)
            @views ρ_buf_up = ρ_buf[1, :]
            @views ρ_buf_down = ρ_buf[2, :]
            optimized_eval_density!(ρ_buf_up, discretization, DUP, X)
            optimized_eval_density!(ρ_buf_down, discretization, DDOWN, X)
            evaluate_vrho!(model; vrho = vρ_buf2, rho = ρ_buf, cache = vρ_buf)
            @views vρ_buf2down = vρ_buf2[2, :]
            Y .= vρ_buf2down
        end

        weightUP = FunWeight(_weight_up!; is_inplace = true, is_vectorized = true)
        weightDOWN = FunWeight(_weight_down!; is_inplace = true, is_vectorized = true)
        fill!(VxcUP, 0)
        fill!(VxcDOWN, 0)
        fill_mass_matrix!(basis, VxcUP; weight = weightUP, method = fem_integration_method)
        fill_mass_matrix!(
            basis, VxcDOWN; weight = weightDOWN, method = fem_integration_method)
        VxcUP .= (VxcUP .+ VxcUP') ./ 2
        VxcDOWN .= (VxcDOWN .+ VxcDOWN') ./ 2
    end
    nothing
end

# ===================================================================
#                       Hamiltonian Matrix
# ===================================================================
"""
    assemble_hamiltonian!(discretization, model, D)

Assemble the Kohn–Sham Hamiltonian matrices from the current density `D`.

This routine updates the internal Hamiltonian blocks stored in `discretization` by
adding the fixed part (Kinetic + Coulomb), the Hartree potential (if enabled), and the
exchange–correlation potential (if present in the model).
"""
function assemble_hamiltonian!(discretization::KSEDiscretization,
                           model::KSEModel,
                           D::AbstractArray{<:Real})
    @unpack H, Hfix, Hartree, VxcUP, VxcDOWN = discretization.ksham

    # Compute hartree matrix
    iszero(model.hartree) || assemble_hartree!(discretization, D, model.N, model.hartree)

    # Compute Exchange-Correlation matrix
    !has_exchcorr(model) || assemble_exc!(discretization, model, D)

    # Build the hamiltonian for each section (l,σ)
    @threads for l in 0:discretization.lₕ
        @views vH = H[:, :, l + 1, 1]
        @. vH = Hfix[l + 1] + VxcUP + Hartree
    end

    if discretization.nspin == 2
        @threads for l in 0:discretization.lₕ
            @views vH = H[:, :, l + 1, 2]
            @. vH = Hfix[l + 1] + VxcDOWN + Hartree
        end
    end
    nothing
end
