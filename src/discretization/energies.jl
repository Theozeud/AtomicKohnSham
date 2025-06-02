#--------------------------------------------------------------------
#                         TOTAL ENERGY
#--------------------------------------------------------------------
function compute_total_energy( discretization::KSEDiscretization, 
                                model::KSEModel,
                                D::AbstractMatrix{<:Real}, 
                                n::AbstractMatrix{<:Real},
                                ϵ::AbstractMatrix{<:Real})
    @unpack Rmax, matrices, exc = discretization
    @unpack VxcUP, VxcDOWN = matrices
    @tensor energy = n[l,n] * ϵ[l,n] 
    if isthereExchangeCorrelation(model)
        energy_exc = compute_exchangecorrelation_energy(discretization, model, D)
        energy_har = compute_hartree_energy(discretization, D)
        if exc == 1
            @tensor energy_correction = Vxc[i,j] * D[i,j]
            return energy - energy_har + energy_exc - energy_correction
        elseif exc == 2
            @views DUP = D[:,:,1]   
            @views DDOWN = D[:,:,2]
            @tensor energy_correctionUP     = VxcUP[i,j] * DUP[i,j]
            @tensor energy_correctionDOWN   = VxcDOWN[i,j] * DDOWN[i,j]
            return energy - energy_har + energy_exc - energy_correctionUP - energy_correctionDOWN
        end
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
    @unpack lₕ, nₕ, exc  = discretization
    elT = eltype(discretization)
    @unpack Kin = discretization.matrices
    energy_kin = zero(elT)
    @inbounds for σ ∈ 1:exc
        @inbounds for l ∈ 1:lₕ+1 
            @inbounds for k ∈ 1:nₕ
                if !iszero(n[l,k,σ])
                    @views Ulkσ = U[:,k,l,σ]
                    energy_kin += n[l,k,σ] * Ulkσ' * Kin[l] * Ulkσ
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
    @unpack lₕ, nₕ, exc  = discretization
    elT = eltype(discretization)
    @unpack Coulomb = discretization.matrices
    energy_cou = zero(elT)
    @inbounds for σ ∈ 1:exc
        @inbounds for l ∈ 1:lₕ+1   
            @inbounds for k ∈ 1:nₕ
                if !iszero(n[l,k,σ])
                    @views Ulkσ = U[:,k,l,σ]
                    energy_cou +=  n[l,k] * Ulkσ' * Coulomb * Ulkσ
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
                                D::AbstractMatrix{<:Real})
    @unpack Rmax, matrices, cache, exc = discretization
    elT = eltype(discretization)
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    if exc == 1
        tensor_matrix_dict!(tmp_B, D, F)
        tmp_C .= A\tmp_B
        @tensor Crho = D[i,j] * M₀[i,j]
        return elT(0.5) * (dot(tmp_B,tmp_C) + Crho^2/Rmax)
    elseif exc == 2
        @views DUP = D[:,:,1]   
        @views DDOWN = D[:,:,2]
        tensor_matrix_dict!(tmp_B, DUP, DDOWN, F)
        tmp_C .= A\tmp_B
        @tensor CrhoUP      = DUP[i,j] * M₀[i,j]
        @tensor CrhoDOWN    = DDOWN[i,j] * M₀[i,j]
        Crho = CrhoUP + CrhoDOWN
        return elT(0.5) * (dot(tmp_B,tmp_C) + Crho^2/Rmax)
    end
end


function compute_hartree_mix_energy(discretization::KSEDiscretization, 
                                    D0::AbstractMatrix{<:Real}, 
                                    D1::AbstractMatrix{<:Real})
    @unpack Rmax, matrices, cache, exc = discretization
    elT = eltype(discretization)
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    if exc == 1
        tensor_matrix_dict!(tmp_B, D0, F)
        tmp_C .= A\tmp_B
        tensor_matrix_dict!(tmp_B, D1, F)
        @tensor Crho0 = D0[i,j] * M₀[i,j]
        @tensor Crho1 = D1[i,j] * M₀[i,j]
        return elT(0.5) * (dot(tmp_B,tmp_C) + Crho0*Crho1/Rmax)
    elseif exc == 2
        @views D0UP     = D0[:,:,1]   
        @views D0DOWN   = D0[:,:,2]
        @views D1UP     = D1[:,:,1]   
        @views D1DOWN   = D1[:,:,2]
        tensor_matrix_dict!(tmp_B, D0UP, D0DOWN, F)
        tmp_C .= A\tmp_B
        tensor_matrix_dict!(tmp_B, D1UP, D1DOWN, F)
        @tensor Crho0UP     = D0UP[i,j] * M₀[i,j]
        @tensor Crho0DOWN   = D0DOWN[i,j] * M₀[i,j]
        @tensor Crho1UP     = D1UP[i,j] * M₀[i,j]
        @tensor Crho1DOWN   = D1DOWN[i,j] * M₀[i,j]
        Crho0 = Crho0UP + Crho0DOWN
        Crho1 = Crho1UP + Crho1DOWN
        return elT(0.5) * (dot(tmp_B,tmp_C) + Crho0*Crho1/Rmax)
    end
end


#--------------------------------------------------------------------
#                  EXCHANGE CORRELATION ENERGY
#--------------------------------------------------------------------
function compute_exchangecorrelation_energy(discretization::KSEDiscretization, 
                                            model::KSEModel, 
                                            D::AbstractMatrix{<:Real})
    @unpack Rmax, exc = discretization
    if exc ==1
        ρ(x) = compute_density(discretization, D, x)
        f1(x,p) = exc(model.exc, ρ(x)) * x^2
        prob = IntegralProblem(f1, (zero(Rmax), Rmax))
        return 4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
    elseif exc == 2
        @views DUP = D[:,:,1]   
        @views DDOWN = D[:,:,2]
        ρUP(x) = compute_densityUP(discretization, DUP, x)
        ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
        f2(x,p) = exc(model.exc, ρUP(x), ρDOWN(x)) * x^2
        prob = IntegralProblem(f2, (zero(Rmax), Rmax))
        return 4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
    end
end


function compute_kinetic_correlation_energy!(discretization::KSEDiscretization, 
                                             model::KSEModel, 
                                             D::AbstractMatrix{<:Real})
    @unpack Rmax, exc = discretization
    elT = eltype(discretization)
    if exc == 1
        return zero(elT)
    end
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
    ρ(x) = ρDOWN(x) + ρUP(x)
    ξ(x) = (ρUP(x) - ρDOWN(x))/ρ(x)
    rs(x) = (3/(4π * ρ(x)))^(1/3)
    tc(x,p) =  -4 * ec(model.exc, ρUP(x), ρDOWN(x)) * ρ(x) * x^2 + 3 * x^2 * ( ρUP(x)* vcUP(model.exc, ρUP(x), ρDOWN(x))+ ρDOWN(x) * vcDOWN(model.exc, ρUP(x), ρDOWN(x))) #
    prob = IntegralProblem(tc, (zero(Rmax),Rmax))
    solver.energy_kincor = 4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
    nothing
end
