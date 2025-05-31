#--------------------------------------------------------------------
#                          Kinetic Matrix
#--------------------------------------------------------------------
function kinetic_matrix!(discretization::LDADiscretization)
    @unpack A, M₋₂, Kin = discretization.matrices
    for l ∈ 0:discretization.lₕ
        #@views vkin = Kin[l+1,:,:]
        @. Kin[l+1] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end


#--------------------------------------------------------------------
#                          Coulomb Matrix
#--------------------------------------------------------------------
function coulomb_matrix!(discretization::LDADiscretization, model::KohnShamExtended)
    @unpack M₋₁, Coulomb = discretization.matrices
    Coulomb .= - model.z .* M₋₁
    nothing
end


#--------------------------------------------------------------------
#                          Hartree Matrix
#--------------------------------------------------------------------
function tensor_matrix_dict!(B::AbstractVector{<:Real}, D::AbstractMatrix{<:Real}, F::Dict{Tuple{Int,Int,Int},<:Real})
    fill!(B, zero(eltype(B)))
    @inbounds for ((i, j, m), F_ijm) ∈ F
        B[m] += D[i, j] * F_ijm
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


function hartree_matrix!(discretization::LDADiscretization, D::AbstractMatrix{<:Real}, coeff::Real = true)
    @unpack Rmax, matrices, cache = discretization
    @unpack A, M₀, F, Hartree = matrices
    @unpack tmp_MV, tmp_B, tmp_C = cache
    tensor_matrix_dict!(tmp_B, D, F)
    tmp_C .= A\tmp_B
    @tensor newCrho = D[i,j] * M₀[i,j]
    tensor_vector_dict!(tmp_MV, tmp_C, F)
    @. Hartree = tmp_MV + newCrho/Rmax * M₀
    @. Hartree .*= coeff
    Hartree .= (Hartree .+ Hartree') ./2
    nothing
end


#--------------------------------------------------------------------
#                   Exchange Correlation Matrix
#--------------------------------------------------------------------
function exchange_corr_matrix!( discretization::LDADiscretization, 
                                model::KohnShamExtended, 
                                D::AbstractMatrix{<:Real})
    @unpack matrices, basis = discretization
    @unpack Vxc = matrices
    ρ(x) = compute_density(discretization, D, x)
    weight = FunWeight(x -> vxc(model.exc, ρ(x)))
    fill!(Vxc, zero(eltype(Vxc))) 
    fill_mass_matrix!(basis, Vxc; weight = weight)
    Vxc .= (Vxc .+ Vxc') ./2
    nothing
end