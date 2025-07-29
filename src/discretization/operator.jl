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


function hartree_matrix!(discretization::KSEDiscretization, 
                         D::AbstractArray{<:Real}, 
                         coeff::Real = true)
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
                                D::AbstractArray{<:Real})
    @unpack matrices, basis, n_spin, fem_integration_method, cache = discretization
    @unpack tmp_ρ, tmp_vρ, tmp_vρ2 = cache
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
        function _weight_up!(Y::AbstractVector, X::AbstractVector)
            @views tmp_ρ_up      = tmp_ρ[1,:]
            @views tmp_ρ_down    = tmp_ρ[2,:]
            optimized_eval_density!(tmp_ρ_up, discretization, DUP, X)
            optimized_eval_density!(tmp_ρ_down, discretization, DDOWN, X)
            evaluate_vrho!(model; vrho=tmp_vρ2, rho=tmp_ρ, cache=tmp_vρ)
            @views tmp_vρ2up = tmp_vρ2[1,:]
            Y .= tmp_vρ2up
        end

        function _weight_down!(Y::AbstractVector, X::AbstractVector)
            @views tmp_ρ_up      = tmp_ρ[1,:]
            @views tmp_ρ_down    = tmp_ρ[2,:]
            optimized_eval_density!(tmp_ρ_up, discretization, DUP, X)
            optimized_eval_density!(tmp_ρ_down, discretization, DDOWN, X)
            evaluate_vrho!(model; vrho=tmp_vρ2, rho=tmp_ρ, cache=tmp_vρ)
            @views tmp_vρ2down = tmp_vρ2[2,:]
            Y .= tmp_vρ2down
        end

        weightUP    = FunWeight(_weight_up!; is_inplace = true, is_vectorized = true)
        weightDOWN  = FunWeight(_weight_down!; is_inplace = true, is_vectorized = true)
        fill!(VxcUP, 0)
        fill!(VxcDOWN, 0)
        fill_mass_matrix!(basis, VxcUP; weight=weightUP, method=fem_integration_method)
        fill_mass_matrix!(basis, VxcDOWN; weight=weightDOWN, method=fem_integration_method)
        VxcUP   .= (VxcUP .+ VxcUP') ./2
        VxcDOWN .= (VxcDOWN .+ VxcDOWN') ./2
    end
    nothing
end
