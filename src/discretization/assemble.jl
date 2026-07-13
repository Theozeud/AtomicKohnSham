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
    @views DIb = D[Ib, Ib]
    @views QgenxIg = Qgenxreshape[Ig, Ig, :]
    @tensor ρ[k] = DIb[i, j] * QgenxIg[i, j, k]
    @.ρ /= X^2
    @. ρ *= 1/4π
    ρ
end

"""
Assemble the exchange–correlation potential operators from the current density `D`.
"""
function assemble_exc!(discretization::KSEDiscretization, model::KSEModel,
                      D::AbstractArray{<:Real})
    @unpack basis, nspin, fem_integration_method, cache = discretization
    @unpack ρ_buf, vρ_buf, vρ_buf2 = cache.excw
    @unpack VxcUP, VxcDOWN = discretization.ksham
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
