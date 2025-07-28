#--------------------------------------------------------------------
#                          Kinetic Matrix
#--------------------------------------------------------------------
function kinetic_matrix!(discretization::KSEDiscretization)
    @unpack A, M₋₂, Kin = discretization.matrices
    for l ∈ 0:discretization.lₕ
        @. Kin[l+1] =  1/2 * (A + l*(l+1)*M₋₂)
    end
    nothing
end


#--------------------------------------------------------------------
#                          Coulomb Matrix
#--------------------------------------------------------------------
function coulomb_matrix!(discretization::KSEDiscretization, model::KSEModel)
    @unpack M₋₁, Coulomb = discretization.matrices
    Coulomb .= - model.z .* M₋₁
    nothing
end


#--------------------------------------------------------------------
#                          Hartree Matrix
#--------------------------------------------------------------------
function tensor_matrix_dict!(B::AbstractVector{<:Real},
                             D::AbstractMatrix{<:Real},
                             F::Dict{Tuple{Int,Int,Int},<:Real})
    fill!(B, zero(eltype(B)))
    @inbounds for ((i, j, m), F_ijm) ∈ F
        B[m] += D[i, j] * F_ijm
    end
    nothing
end


function tensor_matrix_dict!(B::AbstractVector{<:Real},
                             DUP::AbstractMatrix{<:Real},
                             DDOWN::AbstractMatrix{<:Real},
                             F::Dict{Tuple{Int,Int,Int},<:Real})
    fill!(B, zero(eltype(B)))
    @inbounds for ((i, j, m), F_ijm) ∈ F
        B[m] += (DUP[i, j] + DDOWN[i,j]) * F_ijm
    end
    nothing
end


function tensor_vector_dict!(B::AbstractMatrix{<:Real}, D::AbstractVector{<:Real}, F::Dict{Tuple{Int,Int,Int},<:Real})
    fill!(B, zero(eltype(B)))
    @inbounds for ((i, j, m), F_ijm) ∈ F
        B[i,j] += D[m] * F_ijm
    end
    nothing
end


function hartree_matrix!(discretization::KSEDiscretization, D::AbstractMatrix{<:Real}, coeff::Real = true)
    @unpack Rmax, matrices, cache, n_spin = discretization
    @unpack A, M₀, F, Hartree = matrices
    @unpack tmp_MV, tmp_B, tmp_C = cache
    if n_spin == 1
        tensor_matrix_dict!(tmp_B, D, F)
        tmp_C .= A\tmp_B
        @tensor newCrho = D[i,j] * M₀[i,j]
    else
        @views DUP = D[:,:,1]
        @views DDOWN = D[:,:,2]
        tensor_matrix_dict!(tmp_B, DUP, DDOWN, F)
        tmp_C .= A\tmp_B
        @tensor newCrho = DUP[i,j] * M₀[i,j] + DDOWN[i,j] * M₀[i,j]
    end
    tensor_vector_dict!(tmp_MV, tmp_C, F)
    @. Hartree = tmp_MV + newCrho/Rmax * M₀
    @. Hartree .*= coeff
    Hartree .= (Hartree .+ Hartree') ./2
    nothing
end


#--------------------------------------------------------------------
#                   Exchange Correlation Matrix
#--------------------------------------------------------------------
function exchange_corr_matrix!( discretization::KSEDiscretization,
                                model::KSEModel,
                                D::AbstractMatrix{<:Real})
    @unpack matrices, basis, n_spin, fem_integration_method, cache = discretization
    @unpack tmp_ρ, tmp_vρ = cache
    @unpack VxcUP, VxcDOWN = matrices
    if n_spin == 1
        function _weight!(Y::AbstractVector, X::AbstractVector)
            optimized_eval_density!(tmp_ρ, discretization, D, X)
            evaluate_vrho!(model; vrho=Y, rho=tmp_ρ, cache=tmp_vρ)
        end
        weight = FunWeight(_weight!; is_inplace = true, is_vectorized = true)
        fill!(VxcUP, 0)
        fill_mass_matrix!(basis, VxcUP; weight = weight, method = fem_integration_method)
        VxcUP .= (VxcUP .+ VxcUP') ./2
    else
        @views DUP = D[:,:,1]
        @views DDOWN = D[:,:,2]
        ρUP(x) = compute_densityUP(discretization, DUP, x)
        ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
        weightUP    = FunWeight(x -> vxcUP(model.exc, ρUP(x), ρDOWN(x)))
        weightDOWN  = FunWeight(x -> vxcDOWN(model.exc, ρUP(x), ρDOWN(x)))
        fill!(VxcUP, 0)
        fill!(VxcDOWN, 0)
        fill_mass_matrix!(basis, VxcUP; weight=weightUP, method=fem_integration_method)
        fill_mass_matrix!(basis, VxcDOWN; weight=weightDOWN, method=fem_integration_method)
        VxcUP   .= (VxcUP .+ VxcUP') ./2
        VxcDOWN .= (VxcDOWN .+ VxcDOWN') ./2
    end
    nothing
end
