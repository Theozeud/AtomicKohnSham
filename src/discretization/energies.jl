#--------------------------------------------------------------------
#                         TOTAL ENERGY
#--------------------------------------------------------------------
function compute_total_energy( discretization::KSEDiscretization, 
                                model::KSEModel,
                                D::AbstractMatrix{<:Real}, 
                                n::AbstractMatrix{<:Real},
                                ϵ::AbstractMatrix{<:Real})
    @unpack Rmax, matrices = discretization
    @unpack Vxc = matrices
    @tensor energy = n[l,n] * ϵ[l,n] 
    if isthereExchangeCorrelation(model)
        @tensor energy_correction = Vxc[i,j] * D[i,j]
        energy_exc = compute_exchangecorrelation_energy(discretization, model, D)
        energy_har = compute_hartree_energy(discretization, D)
        return energy - energy_har + energy_exc - energy_correction
    else
        energy_har = compute_hartree_energy(discretization, D)
        return energy = energy - energy_har
    end
    nothing
end


#--------------------------------------------------------------------
#                        KINETIC ENERGY
#--------------------------------------------------------------------
function compute_kinetic_energy(discretization::KSEDiscretization, 
                                U::AbstractArray{<:Real}, 
                                n::AbstractArray{<:Real})
    @unpack lₕ, nₕ, elT  = discretization
    @unpack Kin = discretization.matrices
    energy_kin = zero(elT)
    @inbounds for l ∈ 1:lₕ+1 
        @inbounds for k ∈ 1:nₕ
            if !iszero(n[l,k])
                @views Ulk = U[:,k,l]
                energy_kin += n[l,k] * Ulk' * Kin[l] * Ulk
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
    @unpack lₕ, nₕ, elT  = discretization
    @unpack Coulomb = discretization.matrices
    energy_cou = zero(elT)
    @inbounds for l ∈ 1:lₕ+1   
        @inbounds for k ∈ 1:nₕ
            if !iszero(n[l,k])
                @views Ulk = U[:,k,l]
                energy_cou +=  n[l,k] * Ulk' * Coulomb * Ulk
            end
        end
    end
    return energy_cou
end


#--------------------------------------------------------------------
#                        HARTREE ENERGY
#--------------------------------------------------------------------
function compute_hartree_energy(discretization::KSEDiscretization, 
                                D::AbstractMatrix{<:Real})
    @unpack Rmax, elT, matrices, cache = discretization
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    tensor_matrix_dict!(tmp_B, D, F)
    tmp_C .= A\tmp_B
    @tensor Crho = D[i,j] * M₀[i,j]
    return elT(0.5) * (dot(tmp_B,tmp_C) + Crho^2/Rmax)
end


function compute_hartree_mix_energy(discretization::KSEDiscretization, 
                                    D0::AbstractMatrix{<:Real}, 
                                    D1::AbstractMatrix{<:Real})
    @unpack Rmax, elT, matrices, cache = discretization
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    tensor_matrix_dict!(tmp_B, D0, F)
    tmp_C .= A\tmp_B
    tensor_matrix_dict!(tmp_B, D1, F)
    @tensor Crho0 = D0[i,j] * M₀[i,j]
    @tensor Crho1 = D1[i,j] * M₀[i,j]
    return elT(0.5) * (dot(tmp_B,tmp_C) + Crho0*Crho1/Rmax)
end


#--------------------------------------------------------------------
#                  EXCHANGE CORRELATION ENERGY
#--------------------------------------------------------------------
function compute_exchangecorrelation_energy(discretization::KSEDiscretization, 
                                            model::KSEModel, 
                                            D::AbstractMatrix{<:Real})
    @unpack Rmax = discretization
    ρ(x) = compute_density(discretization, D, x)
    f(x,p) = exc(model.exc, ρ(x)) * x^2
    prob = IntegralProblem(f, (zero(Rmax), Rmax))
    4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
end