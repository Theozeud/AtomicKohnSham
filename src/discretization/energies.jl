#--------------------------------------------------------------------
#                         TOTAL ENERGY
#--------------------------------------------------------------------
function compute_total_energy(discretization::KSEDiscretization,
        model::KSEModel,
        D::AbstractArray{<:Real},
        n::AbstractArray{<:Real},
        ϵ::AbstractArray{<:Real})
    @unpack Rmax, matrices, n_spin = discretization
    @unpack VxcUP, VxcDOWN = matrices
    if has_exchcorr(model)
        energy_exc = compute_exchangecorrelation_energy(discretization, model, D)
        energy_har = compute_hartree_energy(discretization, D)
        if n_spin == 1
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
        if n_spin == 1
            @tensor energy = n[l, k] * ϵ[l, k]
            energy_har = compute_hartree_energy(discretization, D)
            return energy = energy - energy_har
        elseif n_spin == 2
            @tensor energy = n[l, k, σ] * ϵ[l, k, σ]
            energy_har = compute_hartree_energy(discretization, D)
            return energy = energy - energy_har
        end
    end
    nothing
end

#--------------------------------------------------------------------
#                        KINETIC ENERGY
#--------------------------------------------------------------------
function compute_kinetic_energy(discretization::KSEDiscretization,
        U::AbstractArray{<:Real},
        n::AbstractArray{<:Real})
    @unpack lₕ, nₕ, n_spin = discretization
    elT = eltype(discretization)
    @unpack Kin = discretization.matrices
    energy_kin = zero(elT)
    @inbounds for σ in 1:n_spin
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

#--------------------------------------------------------------------
#                        COULOMB ENERGY
#--------------------------------------------------------------------
function compute_coulomb_energy(discretization::KSEDiscretization,
        U::AbstractArray{<:Real},
        n::AbstractArray{<:Real})
    @unpack lₕ, nₕ, n_spin = discretization
    elT = eltype(discretization)
    @unpack Coulomb = discretization.matrices
    energy_cou = zero(elT)
    @inbounds for σ in 1:n_spin
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

#--------------------------------------------------------------------
#                        HARTREE ENERGY
#--------------------------------------------------------------------
function compute_hartree_energy(discretization::KSEDiscretization,
        D::AbstractArray{<:Real})
    @unpack Rmax, matrices, cache, n_spin = discretization
    elT = eltype(discretization)
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    if n_spin == 1
        tensor_matrix_dict!(tmp_B, D, F)
        tmp_C .= A\tmp_B
        @tensor Crho = D[i, j] * M₀[i, j]
        return elT(0.5) * (dot(tmp_B, tmp_C) + Crho^2/Rmax)
    else
        @views DUP = D[:, :, 1]
        @views DDOWN = D[:, :, 2]
        tensor_matrix_dict!(tmp_B, DUP, DDOWN, F)
        tmp_C .= A\tmp_B
        @tensor CrhoUP = DUP[i, j] * M₀[i, j]
        @tensor CrhoDOWN = DDOWN[i, j] * M₀[i, j]
        Crho = CrhoUP + CrhoDOWN
        return elT(0.5) * (dot(tmp_B, tmp_C) + Crho^2/Rmax)
    end
end

function compute_hartree_mix_energy(discretization::KSEDiscretization,
        D0::AbstractArray{<:Real},
        D1::AbstractArray{<:Real})
    @unpack Rmax, matrices, cache, n_spin = discretization
    elT = eltype(discretization)
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    if n_spin == 1
        tensor_matrix_dict!(tmp_B, D0, F)
        tmp_C .= A\tmp_B
        tensor_matrix_dict!(tmp_B, D1, F)
        @tensor Crho0 = D0[i, j] * M₀[i, j]
        @tensor Crho1 = D1[i, j] * M₀[i, j]
        return elT(0.5) * (dot(tmp_B, tmp_C) + Crho0*Crho1/Rmax)
    else
        @views D0UP = D0[:, :, 1]
        @views D0DOWN = D0[:, :, 2]
        @views D1UP = D1[:, :, 1]
        @views D1DOWN = D1[:, :, 2]
        tensor_matrix_dict!(tmp_B, D0UP, D0DOWN, F)
        tmp_C .= A\tmp_B
        tensor_matrix_dict!(tmp_B, D1UP, D1DOWN, F)
        @tensor Crho0UP = D0UP[i, j] * M₀[i, j]
        @tensor Crho0DOWN = D0DOWN[i, j] * M₀[i, j]
        @tensor Crho1UP = D1UP[i, j] * M₀[i, j]
        @tensor Crho1DOWN = D1DOWN[i, j] * M₀[i, j]
        Crho0 = Crho0UP + Crho0DOWN
        Crho1 = Crho1UP + Crho1DOWN
        return elT(0.5) * (dot(tmp_B, tmp_C) + Crho0*Crho1/Rmax)
    end
end

#--------------------------------------------------------------------
#                  EXCHANGE CORRELATION ENERGY
#--------------------------------------------------------------------
function compute_exchangecorrelation_energy(discretization::KSEDiscretization,
        model::KSEModel,
        D::AbstractArray{<:Real})
    @unpack Rmax, n_spin, fem_integration_method = discretization
    @unpack fx, fx2, fy, wy, y, shiftx = fem_integration_method
    if n_spin == 1
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

function compute_kinetic_correlation_energy!(discretization::KSEDiscretization,
        model::KSEModel,
        D::AbstractArray{<:Real})
    @unpack Rmax, n_spin = discretization
    elT = eltype(discretization)
    if n_spin == 1
        return zero(elT)
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
