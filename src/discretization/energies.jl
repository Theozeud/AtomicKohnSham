# ===================================================================
#                         TOTAL ENERGY
# ===================================================================
"""
Compute the Kohn–Sham total energy from eigenvalues, occupations, and density.
"""
function compute_total_energy(discretization::KSEDiscretization, model::KSEModel,
                             D::AbstractArray{<:Real}, n::AbstractArray{<:Real},
                             ϵ::AbstractArray{<:Real})
    @unpack VxcUP, VxcDOWN = discretization.ksham
    if has_exchcorr(model)
        energy_exc = compute_exc_energy(discretization, model, D)
        energy_har = compute_hartree_energy(discretization, D)
        if discretization.nspin == 1
            @tensor energy = n[l, k] * ϵ[l, k]
            @tensor energy_correction = VxcUP[i, j] * D[i, j]
            return energy - energy_har + energy_exc - energy_correction
        else
            @tensor energy = n[l, k, σ] * ϵ[l, k, σ]
            @views DUP = D[:, :, 1]
            @views DDOWN = D[:, :, 2]
            @tensor energy_correctionUP = VxcUP[i, j] * DUP[i, j]
            @tensor energy_correctionDOWN = VxcDOWN[i, j] * DDOWN[i, j]
            return energy - energy_har + energy_exc - energy_correctionUP -
                   energy_correctionDOWN
        end
    else
        if discretization.nspin == 1
            @tensor energy = n[l, k] * ϵ[l, k]
            energy_har = compute_hartree_energy(discretization, D)
            return energy = energy - energy_har
        else
            @tensor energy = n[l, k, σ] * ϵ[l, k, σ]
            energy_har = compute_hartree_energy(discretization, D)
            return energy = energy - energy_har
        end
    end
    nothing
end

# ===================================================================
#                        KINETIC ENERGY
# ===================================================================
"""
Compute the kinetic energy from orbital coefficients and occupations.
"""
function compute_kinetic_energy(discretization::KSEDiscretization{T},
                               U::AbstractArray{<:Real}, n::AbstractArray{<:Real}) where T
    @unpack lₕ, nₕ, nspin = discretization
    @unpack Kin = discretization.ksham
    energy_kin = zero(T)
    @inbounds for σ in 1:nspin
        @inbounds for l in 1:(lₕ + 1)
            @inbounds for k in 1:nₕ
                if !iszero(n[l, k, σ])
                    @views Ulkσ = U[:, k, l, σ]
                    energy_kin += n[l, k, σ] * Ulkσ' * Kin[l] * Ulkσ
                end
            end
        end
    end
    return energy_kin
end

# ===================================================================
#                        COULOMB ENERGY
# ===================================================================
"""
Compute the nuclear Coulomb (external potential) energy from orbitals and occupations.
"""
function compute_coulomb_energy(discretization::KSEDiscretization{T},
                               U::AbstractArray{<:Real},n::AbstractArray{<:Real}) where T
    @unpack lₕ, nₕ, nspin = discretization
    @unpack Coulomb = discretization.ksham
    energy_cou = zero(T)
    @inbounds for σ in 1:nspin
        @inbounds for l in 1:(lₕ + 1)
            @inbounds for k in 1:nₕ
                if !iszero(n[l, k, σ])
                    @views Ulkσ = U[:, k, l, σ]
                    energy_cou += n[l, k, σ] * Ulkσ' * Coulomb * Ulkσ
                end
            end
        end
    end
    return energy_cou
end

# ===================================================================
#                        HARTREE ENERGY
# ===================================================================
"""
Compute the Hartree (electrostatic) energy associated with the density matrix `D`.
"""

function compute_hartree_energy(discretization::KSEDiscretization, D::AbstractArray{<:Real})
    @unpack Rmax, cache, nspin, N = discretization
    @unpack A, F = discretization.femops
    @unpack B, W = cache.hartw
    # HartreeWorkspace allocates B/W with length 0 when model.hartree == 0 (see
    # discretization.jl); A\B below would otherwise throw a DimensionMismatch.
    isempty(B) && return zero(eltype(D))
    if nspin == 1
        tensor_matrix_dict!(B, D, F)
        W .= A\B
        return (dot(B, W) + N^2/Rmax)/2
    else
        @views DUP = D[:, :, 1]
        @views DDOWN = D[:, :, 2]
        tensor_matrix_dict!(B, DUP, DDOWN, F)
        W .= A\B
        return (dot(B, W) + N^2/Rmax)/2
    end
end

"""
Compute the mixed Hartree energy D(ρ₀,ρ₁) associated with two densities.
"""
function compute_hartree_mix_energy(discretization::KSEDiscretization,
                                   D0::AbstractArray{<:Real}, D1::AbstractArray{<:Real})
    @unpack Rmax, cache, nspin, N = discretization
    @unpack A, F = discretization.femops
    @unpack B, W = cache.hartw
    # See compute_hartree_energy: B/W are empty when model.hartree == 0.
    isempty(B) && return zero(eltype(D0))
    if nspin == 1
        tensor_matrix_dict!(B, D0, F)
        W .= A\B
        tensor_matrix_dict!(B, D1, F)
        return (dot(B, W) + N^2/Rmax)/2
    else
        @views D0UP = D0[:, :, 1]
        @views D0DOWN = D0[:, :, 2]
        @views D1UP = D1[:, :, 1]
        @views D1DOWN = D1[:, :, 2]
        tensor_matrix_dict!(B, D0UP, D0DOWN, F)
        W .= A\B
        tensor_matrix_dict!(B, D1UP, D1DOWN, F)
        return (dot(B, W) + N^2/Rmax)/2
    end
end

# ===================================================================
#                  EXCHANGE CORRELATION ENERGY
# ===================================================================
"""
Compute the exchange–correlation energy for the density matrix `D`.
"""

function compute_exc_energy(discretization::KSEDiscretization{T}, model::KSEModel,
                           D::AbstractArray{<:Real}) where T
    @unpack Rmax, nspin, fem_integration_method, cache, basis = discretization
    @unpack fx, fx2, fy, wy, y, shiftx = fem_integration_method
    if !has_exchcorr(model)
        return zero(T)
    end
    if has_gga(model)
        # Unlike LDA (whose ρ-only integrand is C⁰ and integrates fine on the
        # single global Gauss-Legendre grid `y`), GGA's σ=ρ'² integrand is
        # genuinely discontinuous at every mesh node (this basis is only C⁰,
        # so φᵢ' jumps across cell boundaries) -- a single high-order rule
        # spanning all cells converges very poorly across those jumps. So,
        # like the matrix assembly (fill_mass_matrix!/fill_mixed_mass_matrix!),
        # integrate cell-by-cell using fem_integration_method's per-cell
        # reference quadrature (x,w), which never asks one quadrature panel to
        # span a discontinuity. (Confirmed empirically: at fixed global
        # npoints, an analytic-vs-FD check of assemble_exc!'s VxcUP against
        # this per-cell compute_exc_energy converges properly, whereas the
        # single-global-grid version required a ~60x finer grid to agree.)
        @unpack dρ_buf, σ_buf = cache.excw
        @unpack x, w = fem_integration_method
        energy = zero(T)
        if nspin == 1
            for k in cellrange(basis.mesh)
                eldata = getelement(basis, k, :M)
                @. shiftx = eldata.invϕ[1] * x + eldata.invϕ[2]
                optimized_eval_density!(fx, discretization, D, shiftx)
                optimized_eval_density_gradient!(dρ_buf, discretization, D, fx, shiftx)
                @. σ_buf = dρ_buf^2
                evaluate_zk!(model; rho = fx, sigma = σ_buf, zk = fy, cache = dρ_buf)
                @. fx = fx * fy * shiftx^2
                energy += eldata.invϕ[1] * dot(w, fx)
            end
            return 4π * energy
        else
            @views ρup = fx2[1, :]
            @views ρdown = fx2[2, :]
            @views dρup = dρ_buf[1, :]
            @views dρdown = dρ_buf[2, :]
            @views DUP = D[:, :, 1]
            @views DDOWN = D[:, :, 2]
            @views σuu = σ_buf[1, :]
            @views σud = σ_buf[2, :]
            @views σdd = σ_buf[3, :]
            for k in cellrange(basis.mesh)
                eldata = getelement(basis, k, :M)
                @. shiftx = eldata.invϕ[1] * x + eldata.invϕ[2]
                optimized_eval_density!(ρup, discretization, DUP, shiftx)
                optimized_eval_density!(ρdown, discretization, DDOWN, shiftx)
                optimized_eval_density_gradient!(dρup, discretization, DUP, ρup, shiftx)
                optimized_eval_density_gradient!(dρdown, discretization, DDOWN, ρdown, shiftx)
                @. σuu = dρup^2
                @. σud = dρup * dρdown
                @. σdd = dρdown^2
                # cache must be a genuine Array, not a view (dρup/dρdown are
                # SubArrays of dρ_buf's rows): Libxc's low-level evaluate!
                # methods require zk/vrho/... ::Union{Array{Float64},Ptr}
                # exactly. `fx` is unused until the line below, so reuse it.
                evaluate_zk!(model; rho = fx2, sigma = σ_buf, zk = fy, cache = fx)
                @. fx = ρup + ρdown
                @. fx = fx * fy * shiftx^2
                energy += eldata.invϕ[1] * dot(w, fx)
            end
            return 4π * energy
        end
    end
    if nspin == 1
        eval_density!(fx, discretization, D, y)
        evaluate_zk!(model; rho = fx, zk = fy, cache = shiftx)
        @. fx = fx * fy * y^2
        return 4π * dot(wy, fx)
    else
        @views ρup = fx2[1, :]
        @views ρdown = fx2[2, :]
        @views DUP = D[:, :, 1]
        @views DDOWN = D[:, :, 2]
        eval_density!(ρup, discretization, DUP, y)
        eval_density!(ρdown, discretization, DDOWN, y)
        evaluate_zk!(model; rho = fx2, zk = fy, cache = shiftx)
        @. fx = ρup + ρdown
        @. fx = fx * fy * y^2
        return 4π * dot(wy, fx)
    end
end

#=
function compute_kinetic_correlation_energy!(discretization::KSEDiscretization,
        model::KSEModel,
        D::AbstractArray{<:Real})
    @unpack Rmax, nspin = discretization
    if nspin == 1
        return zero(Rmax)
    end
    @views DUP = D[:, :, 1]
    @views DDOWN = D[:, :, 2]
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
    ρ(x) = ρDOWN(x) + ρUP(x)
    ξ(x) = (ρUP(x) - ρDOWN(x))/ρ(x)
    rs(x) = (3/(4π * ρ(x)))^(1/3)
    tc(x,
        p) = -4 * ec(model.exc, ρUP(x), ρDOWN(x)) * ρ(x) * x^2 +
             3 * x^2 *
             (ρUP(x) *
              vcUP(model.exc, ρUP(x),
                 ρDOWN(x))+ρDOWN(x) * vcDOWN(model.exc, ρUP(x), ρDOWN(x)))
    prob = IntegralProblem(tc, (zero(Rmax), Rmax))
    solver.energy_kincor = 4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
    nothing
end
=#
