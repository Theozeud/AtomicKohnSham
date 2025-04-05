#####################################################################
#                          LSDA Cache
#####################################################################

mutable struct LSDAMatrices{T<:Real, 
                            typeMatrix <: AbstractMatrix{<:Real}, 
                            typeVectorMatrix <: Vector{<:AbstractMatrix{<:Real}}}
    # FEM MATRICES
    A::typeMatrix                   # Matrix of Qᵢ'Qⱼ'
    M₀::Matrix{T}                   # Matrix of QᵢQⱼ
    M₋₁::typeMatrix                 # Matrix of 1/x QᵢQⱼ
    M₋₂::typeMatrix                 # Matrix of 1/x² QᵢQⱼ
    F::Dict{Tuple{Int,Int,Int},T}   # Tensor of 1/x QᵢQⱼQₖ            
    # MATRICES COMPOSING THE HAMILTONIAN
    H::Array{T,4}                   # Hamiltonian
    Kin::typeVectorMatrix           # Kinetic Matrix
    Coulomb::typeMatrix             # Coulomb Matrix
    Hfix::typeVectorMatrix          # (Kinetic + Coulomb) Matrix        
    Hartree::typeMatrix             # Hartree Matrix 
    VxcUP::typeMatrix               # Echange-correlation Matrix UP
    VxcDOWN::typeMatrix             # Echange-correlation Matrix DOWN
end

mutable struct LSDACache{T <: Real}
    tmp_MV::Matrix{T}               # Store the contraction F:C
    tmp_B::Vector{T}                # Matrix of 4πρQᵢ  
    tmp_C::Vector{T}                # Matrix for V(ρ) - solution of the Gauss Electrostatic law
    tmp_vect::Vector{T}             # Vector to store temprary vector
end

function create_cache_lsda(lₕ::Int, Nₕ::Int, T::Type)
    # FEM MATRICES
    A           = spzeros(T, Nₕ, Nₕ) 
    M₀          = zeros(T, Nₕ, Nₕ)
    M₋₁         = spzeros(T, Nₕ, Nₕ)
    M₋₂         = spzeros(T, Nₕ, Nₕ)
    F           = Dict{Tuple{Int,Int,Int},T}()
    # MATRICES COMPOSING THE HAMILTONIAN
    H           = zeros(T, Nₕ, Nₕ, lₕ+1, 2)
    Kin         = _spzeros(T, Nₕ, Nₕ, lₕ+1)
    Coulomb     = spzeros(T, Nₕ, Nₕ)
    Hfix        = _spzeros(T, Nₕ, Nₕ, lₕ+1)
    Hartree     = spzeros(T, Nₕ, Nₕ)
    VxcUP       = spzeros(T, Nₕ, Nₕ)
    VxcDOWN     = spzeros(T, Nₕ, Nₕ)

    # Initialization of array for temporary stockage of computations 
    tmp_MV          = zeros(T, Nₕ, Nₕ)  
    tmp_B           = zeros(T, Nₕ)
    tmp_C           = zeros(T, Nₕ)
    tmp_vect        = zeros(T, Nₕ)

    LSDAMatrices{T, typeof(A), typeof(Hfix)}(A, M₀, M₋₁, M₋₂, F, H, Kin, Coulomb, Hfix, Hartree, VxcUP, VxcDOWN),  
    LSDACache{T}(tmp_MV, tmp_B, tmp_C, tmp_vect)
end


#####################################################################
#                          LSDA Discretization
#####################################################################


struct LSDADiscretization{T <: Real, T2 <: Real, typeBasis <: Basis} <: KohnShamDiscretization
    lₕ::Int
    nₕ::Int
    Nₕ::Int
    basis::typeBasis
    mesh::Mesh{T}
    Rmin::T
    Rmax::T
    elT::Type
    matrices::LSDAMatrices{T2}
    cache::LSDACache{T2}
    function LSDADiscretization(lₕ::Int, basis::Basis, mesh::Mesh, nₕ::Int = length(basis))
        elT = eltype(basis)
        Nₕ = length(basis)
        new{eltype(mesh), elT, typeof(basis)}(lₕ, nₕ, Nₕ, basis, mesh, first(mesh), last(mesh), elT, create_cache_lsda(lₕ, Nₕ, elT)...)
    end
end

Base.eltype(discretization::LSDADiscretization) = discretization.elT
dim(discretization::LSDADiscretization) = discretization.Nₕ * (discretization.lₕ + 1)

#####################################################################
#                          Init Cache
#####################################################################

function init_cache!(discretization::LSDADiscretization, model::AbstractDFTModel, hartree::Real, integration_method::IntegrationMethod)

    @unpack lₕ, basis, matrices  = discretization
    @unpack A, M₀, M₋₁, M₋₂, F, Kin, Coulomb, Hfix = matrices

    # CREATION OF FEM MATRICES
    fill_stiffness_matrix!(basis, A; method = integration_method)
    fill_mass_matrix!(basis, M₀; method = integration_method)
    fill_mass_matrix!(basis, -1, M₋₁; method = integration_method)
    lₕ == 0 || fill_mass_matrix!(basis, -2, M₋₂; method = integration_method)
    iszero(hartree) || fill_mass_tensor!(basis, -1, F; method = integration_method)

    # CREATION OF THE FIX PART OF THE HAMILTONIAN 
    kinetic_matrix!(discretization)
    coulomb_matrix!(discretization, model)
    for l ∈ 1:lₕ+1
        @. Hfix[l] = Kin[l] + Coulomb
    end
    nothing
end

#####################################################################
#                          Initialization
#####################################################################

init_density(kd::LSDADiscretization)                 = fill(one(kd.elT), kd.Nₕ, kd.Nₕ, 2)  
init_orbitals(kd::LSDADiscretization)                = zeros(kd.elT, kd.Nₕ, kd.nₕ, kd.lₕ+1, 2)
init_orbitals_energy(kd::LSDADiscretization)         = zeros(kd.elT, kd.lₕ+1, kd.nₕ, 2)
init_occupation_number(kd::LSDADiscretization)       = zeros(kd.elT, kd.lₕ+1, kd.nₕ, 2)
init_density_matrix(kd::LSDADiscretization)          = BlockDiagonal([zeros(kd.elT, kd.Nₕ, kd.Nₕ) for i ∈ 1:kd.lₕ+1])

function init_energies(kd::LSDADiscretization, model::KohnShamExtended)
    @unpack elT = kd
    d = Dict(   :Etot => zero(elT),                                     # Total energy 
                :Ekin => zero(elT),                                     # Kinetic energy
                :Ecou => zero(elT),                                     # Coulomb energy
                :Ehar => zero(elT))                                     # Hartree energy     
    !(isthereExchangeCorrelation(model)) ||  (d[:Eexc] = zero(elT))     # Exchange-correlation energy
    !(isthereExchangeCorrelation(model)) ||  (d[:Ekincorr] = zero(elT)) # Kinetic-correlation energy
    d
end

#####################################################################
#               Find Orbital : Solve the eigen problems
#####################################################################

function prepare_eigenvalue_problem!(   discretization::LSDADiscretization, 
                                        model::KohnShamExtended, 
                                        D::AbstractArray{<:Real}, 
                                        hartree::Real = true)

    @unpack H, Hfix, Hartree, VxcUP, VxcDOWN = discretization.matrices

    # COMPUTE HARTREE MATRIX
    iszero(hartree) || hartree_matrix!(discretization, D, hartree)

    # COMPUTE EXCHANGE CORRELATION MATRIX
    !isthereExchangeCorrelation(model) || 
                        exchange_corr_matrix!(discretization, model, D)

    # BUILD THE HAMILTONIAN OF THE lᵗʰ SECTION FOR EACH SPIN σ ∈ {↑,↓}
    @threads for l ∈ 0:discretization.lₕ
        @views vH = H[:,:,l+1,1]
        @. vH = Hfix[l+1] + VxcUP + Hartree
    end
    @threads for l ∈ 0:discretization.lₕ
        @views vH = H[:,:,l+1,2]
        @. vH = Hfix[l+1] + VxcDOWN + Hartree
    end
    nothing
end

function find_orbital!( discretization::LSDADiscretization, 
                        U::AbstractArray{<:Real}, 
                        ϵ::AbstractArray{<:Real})

    @unpack lₕ, nₕ, matrices = discretization
    @unpack M₀, H = matrices

    # SOLVE THE GENERALIZED EIGENVALUE PROBLEM FOR EACH SECTION l
    for σ ∈ 1:2
        for l ∈ 0:lₕ
            @views vH = H[:,:,l+1,σ]
            ϵ[l+1,:,σ], U[:,:,l+1,σ] = solve_generalized_eigenvalue_problem(vH, M₀, nₕ)        
        end
    end
end

#####################################################################
#               Normamisation of eigenvector
#####################################################################

function normalization!(discretization::LSDADiscretization, 
                        U::AbstractArray{<:Real}, 
                        n::AbstractArray{<:Real})
    @unpack M₀ = discretization.matrices
    @unpack lₕ, nₕ = discretization
    @inbounds for k ∈ 1:nₕ
        @inbounds for σ ∈ 1:2
            @inbounds for l ∈ 1:lₕ+1   
                if !iszero(n[l,k,σ])
                    @views Ulk = U[:,k,l,σ]
                    coeff = sqrt(Ulk'*M₀*Ulk)
                    Ulk .= Ulk .* 1.0/coeff
                end
            end
        end
    end
    nothing
end

function normalization!(discretization::LSDADiscretization, U::AbstractArray{<:Real}, l::Int, k::Int, σ::Int)
    @unpack M₀ = discretization.matrices
    @views Ulk = U[:,k,l+1,σ]
    coeff = sqrt(Ulk'*M₀*Ulk)
    Ulk .= Ulk .* 1.0/coeff
    nothing
end

#####################################################################
#                          Kinetic Matrix
#####################################################################

function kinetic_matrix!(discretization::LSDADiscretization)
    @unpack A, M₋₂, Kin = discretization.matrices
    for l ∈ 0:discretization.lₕ
        @. Kin[l+1] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end

#####################################################################
#                          Coulomb Matrix
#####################################################################

function coulomb_matrix!(discretization::LSDADiscretization, model::KohnShamExtended)
    @unpack M₋₁, Coulomb = discretization.matrices
    Coulomb .= - model.z .* M₋₁
    nothing
end

#####################################################################
#                          Hartree Matrix
#####################################################################

function tensor_matrix_dict!(B::AbstractVector{<:Real}, DUP::AbstractMatrix{<:Real}, DDOWN::AbstractMatrix{<:Real}, F::Dict{Tuple{Int,Int,Int},<:Real})
    fill!(B, zero(eltype(B)))
    @inbounds for ((i, j, m), F_ijm) ∈ F
        B[m] += (DUP[i, j] + DDOWN[i,j]) * F_ijm
    end
    nothing
end

function hartree_matrix!(discretization::LSDADiscretization, D::AbstractArray{<:Real}, coeff::Real = true)
    @unpack Rmax, matrices, cache = discretization
    @unpack A, M₀, F, Hartree = matrices
    @unpack tmp_MV, tmp_B, tmp_C = cache
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    tensor_matrix_dict!(tmp_B, DUP, DDOWN, F)
    tmp_C .= A\tmp_B
    @tensor newCrho = DUP[i,j] * M₀[i,j] + DDOWN[i,j] * M₀[i,j]
    tensor_vector_dict!(tmp_MV, tmp_C, F)
    @. Hartree = tmp_MV + newCrho/Rmax * M₀
    @. Hartree .*= coeff
    Hartree .= (Hartree .+ Hartree') ./2
    nothing
end

#####################################################################
#                   Exchange Correlation Matrix
#####################################################################

function exchange_corr_matrix!( discretization::LSDADiscretization, 
                                model::KohnShamExtended, 
                                D::AbstractArray{<:Real})
    @unpack matrices, basis = discretization
    @unpack VxcUP, VxcDOWN = matrices
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
    weightUP    = FunWeight(x -> vxcUP(model.exc, ρUP(x), ρDOWN(x)))
    weightDOWN  = FunWeight(x -> vxcDOWN(model.exc, ρUP(x), ρDOWN(x)))
    fill!(VxcUP, zero(eltype(VxcUP)))
    fill!(VxcDOWN, zero(eltype(VxcDOWN))) 
    fill_mass_matrix!(basis, VxcUP; weight = weightUP)
    fill_mass_matrix!(basis, VxcDOWN; weight = weightDOWN)
    VxcUP   .= (VxcUP .+ VxcUP') ./2
    VxcDOWN .= (VxcDOWN .+ VxcDOWN') ./2
    nothing
end

#####################################################################
#                         TOTAL ENERGY
#####################################################################

function compute_total_energy(  discretization::LSDADiscretization, 
                                model::KohnShamExtended,
                                D::AbstractArray{<:Real}, 
                                n::AbstractArray{<:Real},
                                ϵ::AbstractArray{<:Real})
    @unpack Rmax, matrices = discretization
    @unpack VxcUP, VxcDOWN = matrices
    @tensor energy = n[l,n,σ] * ϵ[l,n,σ] 
    if isthereExchangeCorrelation(model)
        @views DUP = D[:,:,1]   
        @views DDOWN = D[:,:,2]
        @tensor energy_correctionUP     = VxcUP[i,j] * DUP[i,j]
        @tensor energy_correctionDOWN   = VxcDOWN[i,j] * DDOWN[i,j]
        energy_exc = compute_exchangecorrelation_energy(discretization, model, D)
        energy_har = compute_hartree_energy(discretization, D)
        return energy - energy_har + energy_exc - energy_correctionUP - energy_correctionDOWN
    else
        energy_har = compute_hartree_energy(discretization, D)
        return energy = energy - energy_har
    end
    nothing
end

#####################################################################
#                        KINETIC ENERGY
#####################################################################

function compute_kinetic_energy(discretization::LSDADiscretization, 
                                U::AbstractArray{<:Real}, 
                                n::AbstractArray{<:Real})
    @unpack lₕ, nₕ, elT  = discretization
    @unpack Kin = discretization.matrices
    energy_kin = zero(elT)
    @inbounds for σ ∈ 1:2
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

#####################################################################
#                        COULOMB ENERGY
#####################################################################

function compute_coulomb_energy(discretization::LSDADiscretization, 
                                U::AbstractArray{<:Real}, 
                                n::AbstractArray{<:Real})
    @unpack lₕ, nₕ, elT  = discretization
    @unpack Coulomb = discretization.matrices
    energy_cou = zero(elT)
    @inbounds for σ ∈ 1:2
        @inbounds for l ∈ 1:lₕ+1   
            @inbounds for k ∈ 1:nₕ
                if !iszero(n[l,k,σ])
                    @views Ulkσ = U[:,k,l,σ]
                    energy_cou +=  n[l,k,σ] * Ulkσ' * Coulomb * Ulkσ
                end
            end
        end
    end
    return energy_cou
end

#####################################################################
#                        HARTREE ENERGY
#####################################################################

function compute_hartree_energy(discretization::LSDADiscretization, 
                                D::AbstractArray{<:Real})
    @unpack Rmax, elT, matrices, cache = discretization
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    tensor_matrix_dict!(tmp_B, DUP, DDOWN, F)
    tmp_C .= A\tmp_B
    @tensor CrhoUP      = DUP[i,j] * M₀[i,j]
    @tensor CrhoDOWN    = DDOWN[i,j] * M₀[i,j]
    Crho = CrhoUP + CrhoDOWN
    return elT(0.5) * (dot(tmp_B,tmp_C) + Crho^2/Rmax)
end


function compute_hartree_mix_energy(discretization::LSDADiscretization, 
                                    D0::AbstractArray{<:Real}, 
                                    D1::AbstractArray{<:Real})
    @unpack Rmax, elT, matrices, cache = discretization
    @unpack A, F, M₀ = matrices
    @unpack tmp_B, tmp_C = cache
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

#####################################################################
#                  EXCHANGE CORRELATION ENERGY
#####################################################################

function compute_exchangecorrelation_energy(discretization::LSDADiscretization, 
                                            model::KohnShamExtended, 
                                            D::AbstractArray{<:Real})
    @unpack Rmax = discretization
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
    f(x,p) = exc(model.exc, ρUP(x), ρDOWN(x)) * x^2
    prob = IntegralProblem(f, (zero(Rmax), Rmax))
    4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
end

function compute_kinetic_correlation_energy!(   discretization::LSDADiscretization, 
                                                model::KohnShamExtended,
                                                D::AbstractArray{<:Real})
    @unpack Rmax = discretization
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


#####################################################################
#                             Density
#####################################################################

function density!(  discretization::LSDADiscretization, 
                    U::AbstractArray{<:Real}, 
                    n::AbstractArray{<:Real}, 
                    D::AbstractArray{<:Real})
    @unpack lₕ, nₕ, Nₕ, elT  = discretization
    fill!(D, zero(elT))
    @inbounds for k ∈ 1:nₕ
        @inbounds for σ ∈ 1:2
            @inbounds for l ∈ 1:lₕ+1   
                if !iszero(n[l,k,σ])
                    @inbounds for i ∈ 1:Nₕ
                        val = n[l,k,σ] * U[i,k,l,σ] 
                        @inbounds @simd for j ∈ 1:i
                            D[i,j,σ] += val * U[j,k,l,σ]
                        end
                    end
                end
            end
        end
    end
    @inbounds for i in 1:Nₕ
        @inbounds @simd for j in 1:i-1
            D[j,i,:] .= D[i,j,:]
        end
    end
    nothing
end

function compute_density(discretization::LSDADiscretization, D::AbstractArray{<:Real}, x::Real)
    @unpack basis, cache = discretization
    @unpack tmp_vect, tmp_C = cache
    localisation_x = findindex(basis.mesh, x)
    I = basis.cells_to_indices[localisation_x]
    @views eval_basis = tmp_C[I]
    @inbounds for (n,i) ∈ enumerate(I)
        eval_basis[n] = basis(i,x)
    end
    @views tv = tmp_vect[I]
    @views DUPview      = D[I,I,1]
    @views DDOWNview    = D[I,I,2]
    mul!(tv,DUPview .+ DDOWNview,eval_basis)
    return 1/(4π*x^2) * dot(eval_basis,tv)
end

function compute_densityUP(discretization::LSDADiscretization, DUP::AbstractArray{<:Real}, x::Real)
    @unpack basis, cache = discretization
    @unpack tmp_vect, tmp_C = cache
    localisation_x = findindex(basis.mesh, x)
    I = basis.cells_to_indices[localisation_x]
    @views eval_basis = tmp_C[I]
    @inbounds for (n,i) ∈ enumerate(I)
        eval_basis[n] = basis(i,x)
    end
    @views tv = tmp_vect[I]
    @views DUPview = DUP[I,I,1]
    mul!(tv,DUPview,eval_basis)
    return 1/(4π*x^2) * dot(eval_basis,tv)
end

function compute_densityDOWN(discretization::LSDADiscretization, DDOWN::AbstractArray{<:Real}, x::Real)
    @unpack basis, cache = discretization
    @unpack tmp_vect, tmp_C = cache
    localisation_x = findindex(basis.mesh, x)
    I = basis.cells_to_indices[localisation_x]
    @views eval_basis = tmp_C[I]
    @inbounds for (n,i) ∈ enumerate(I)
        eval_basis[n] = basis(i,x)
    end
    @views tv = tmp_vect[I]
    @views DDOWNview = DDOWN[I,I,1]
    mul!(tv,DDOWNview,eval_basis)
    return 1/(4π*x^2) * dot(eval_basis,tv)
end