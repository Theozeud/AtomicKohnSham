#--------------------------------------------------------------------
#                          KSE MATRICES
#--------------------------------------------------------------------
struct KSEMatrices{ T<:Real, 
                    typeMatrix <: AbstractMatrix{<:Real}, 
                    typeVectorMatrix <: Vector{<:AbstractMatrix{<:Real}},
                    typeHam <: AbstractArray{<:Real}}
    # FEM MATRICES
    A::typeMatrix                   # Matrix of Qᵢ'Qⱼ'
    M₀::Matrix{T}                   # Matrix of QᵢQⱼ
    M₋₁::typeMatrix                 # Matrix of 1/x QᵢQⱼ
    M₋₂::typeMatrix                 # Matrix of 1/x² QᵢQⱼ
    F::Dict{Tuple{Int,Int,Int},T}   # Tensor of 1/x QᵢQⱼQₖ 
    S::Matrix{T}                    # √(M₀⁻¹)          
    # MATRICES COMPOSING THE HAMILTONIAN
    H::typeHam                      # Hamiltonian for each l and each spin
    Kin::typeVectorMatrix           # Kinetic Matrix for each l
    Coulomb::typeMatrix             # Coulomb Matrix
    Hfix::typeVectorMatrix          # (Kinetic + Coulomb) Matrix        
    Hartree::typeMatrix             # Hartree Matrix 
    VxcUP::typeMatrix               # Echange-correlation Matrix UP
    VxcDOWN::typeMatrix             # Echange-correlation Matrix DOWN
                                    # For LDA, only VxcUP is used.
end


#--------------------------------------------------------------------
#                           KSE CACHE
#--------------------------------------------------------------------
struct DiscretizationCache{T <: Real}
    tmp_MV::Matrix{T}               # Store the contraction F:C
    tmp_B::Vector{T}                # Matrix of 4πρQᵢ  
    tmp_C::Vector{T}                # Matrix for V(ρ) - solution of the Gauss Electrostatic law
    tmp_vect::Vector{T}             # Vector to store temprary vector
end


function create_cache_discretization(lₕ::Int, Nₕ::Int, T::Type, exc::Int)
    # FEM MATRICES
    A           = spzeros(T, Nₕ, Nₕ) 
    M₀          = zeros(T, Nₕ, Nₕ)
    M₋₁         = spzeros(T, Nₕ, Nₕ)
    M₋₂         = spzeros(T, Nₕ, Nₕ)
    F           = Dict{Tuple{Int,Int,Int},T}()
    S           = zeros(T, Nₕ, Nₕ)
    # MATRICES COMPOSING THE HAMILTONIAN
    H           = flexible_zeros(T, (Nₕ, Nₕ, lₕ+1), exc)
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

    KSEMatrices{T, typeof(A), typeof(Hfix), typeof(H)}(A, M₀, M₋₁, M₋₂, F, S, 
                                            H, Kin, Coulomb, Hfix, Hartree, VxcUP, VxcDOWN),  
    DiscretizationCache{T}(tmp_MV, tmp_B, tmp_C, tmp_vect)
end


#--------------------------------------------------------------------
#                          KSE Discretization
#--------------------------------------------------------------------
"""
    Structure holding the discretization parameters, the discretized operators and the caches variables.

    The discretization of an element of H^1(R^3) is the following :
    u(r,θ,φ) = ∑(0≤l≤→lₕ)∑(-lₕ≤m≤lₕ)∑(1≤n≤Nₕ) unlm Qn(r)/r Yₗᵐ(θ,φ)

    where (Qn) is the FEM basis discretizing the radial part in H^10(R^3),
    (Yₗᵐ) the spherical harmonics and unlm the coefficients in this discretization.

"""
struct KSEDiscretization{T <: Real, 
                         B <: FEMBasis,
                         M <: Mesh} 
    lₕ::Int
    nₕ::Int
    Nₕ::Int
    basis::B
    mesh::M
    Rmin::T
    Rmax::T
    exc::Int
    matrices::KSEMatrices{T}
    cache::DiscretizationCache{T}
    function KSEDiscretization(lₕ::Int, basis::FEMBasis, mesh::Mesh, exc::Symbol, nₕ::Int = length(basis))
        elT = eltype(basis)
        Nₕ = length(basis)
        intexc =  exc == :lda ? 1 : 2
        new{elT, typeof(basis), typeof(mesh)}(lₕ, nₕ, Nₕ, basis, mesh, elT(first(mesh)), elT(last(mesh)),intexc,
                                                create_cache_discretization(lₕ, Nₕ, elT, intexc)...)
    end
end


function init_cache!(discretization::KSEDiscretization, model::KSEModel, integration_method::IntegrationMethod)

    @unpack lₕ, basis, matrices  = discretization
    @unpack A, M₀, M₋₁, M₋₂, F, S, Kin, Coulomb, Hfix = matrices

    # CREATION OF FEM MATRICES
    fill_stiffness_matrix!(basis, A; method = integration_method)
    fill_mass_matrix!(basis, M₀; method = integration_method)
    fill_mass_matrix!(basis, -1, M₋₁; method = integration_method)
    lₕ == 0 || fill_mass_matrix!(basis, -2, M₋₂; method = integration_method)
    iszero(model.hartree) || fill_mass_tensor!(basis, -1, F; method = integration_method)
    S .= sqrt(inv(Symmetric(M₀)))

    # CREATION OF THE FIX PART OF THE HAMILTONIAN 
    kinetic_matrix!(discretization)
    coulomb_matrix!(discretization, model)
    for l ∈ 1:lₕ+1
        @. Hfix[l] = Kin[l] + Coulomb
    end
    nothing
end

#--------------------------------------------------------------------
#                               API
#--------------------------------------------------------------------

Base.eltype(::KSEDiscretization{T}) where T = T
dim(discretization::KSEDiscretization) = discretization.Nₕ * (discretization.lₕ + 1)
multiplicity(::KSEDiscretization) = 2


zero_density(kd::KSEDiscretization)                 = flexible_zeros(eltype(kd), (kd.Nₕ, kd.Nₕ), kd.exc)  
zero_orbitals(kd::KSEDiscretization)                = flexible_zeros(eltype(kd), (kd.Nₕ, kd.nₕ, kd.lₕ+1), kd.exc)
zero_orbitals_energy(kd::KSEDiscretization)         = flexible_zeros(eltype(kd), (kd.lₕ+1, kd.nₕ), kd.exc)
zero_occupation_number(kd::KSEDiscretization)       = flexible_zeros(eltype(kd), (kd.lₕ+1, kd.nₕ), kd.exc)
zero_density_matrix(kd::KSEDiscretization)          = flexible_zeros(eltype(kd), (kd.Nₕ, kd.Nₕ, kd.lₕ+1), kd.exc)

function zero_energies(kd::KSEDiscretization, model::KSEModel)
    elT = eltype(kd)
    d = Dict(   :Etot => zero(elT),                                     # Total energy 
                :Ekin => zero(elT),                                     # Kinetic energy
                :Ecou => zero(elT),                                     # Coulomb energy
                :Ehar => zero(elT))                                     # Hartree energy
    if isthereExchangeCorrelation(model)
        (d[:Eexc] = zero(elT))                                          # Exchange-correlation energy
        kd.exc != :lsda || (d[:Ekincorr] = zero(elT))                   # Kinetic-correlation energy
    end
    d
end

zero_operator(kd::KSEDiscretization)        = zeros(eltype(kd), kd.Nₕ, kd.Nₕ, kd.lₕ+1)
zero_single_operator(kd::KSEDiscretization) = zeros(eltype(kd), kd.Nₕ, kd.Nₕ)


#--------------------------------------------------------------------
#               Find Orbital : Solve the eigen problems
#--------------------------------------------------------------------
function build_hamiltonian!(discretization::KSEDiscretization, 
                            model::KSEModel, 
                            D::AbstractMatrix{<:Real}, 
                            hartree::Real = true)

    @unpack H, Hfix, Hartree, VxcUP = discretization.matrices

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

    if discretization.exc == 2
        @threads for l ∈ 0:discretization.lₕ
            @views vH = H[:,:,l+1,2]
            @. vH = Hfix[l+1] + VxcDOWN + Hartree
        end
    end
    nothing
end


function find_orbital!( discretization::KSEDiscretization, 
                        U::AbstractArray{<:Real}, 
                        ϵ::AbstractMatrix{<:Real})

    @unpack lₕ, nₕ, matrices, exc = discretization
    @unpack S, H = matrices

    # SOLVE THE GENERALIZED EIGENVALUE PROBLEM FOR EACH SECTION l AND FOR EACH σ ∈ {↑,↓}
    for σ ∈ 1:exc
        for l ∈ 0:lₕ
            @views vH = H[:,:,l+1,σ]
            λ, V = eigen(Symmetric(S*vH*S))
            @views ϵ[l+1,:,σ] =  λ[1:nₕ]
            @views U[:,:,l+1,σ] = S*V[:,1:nₕ]        
        end
    end
end


#--------------------------------------------------------------------
#               Normamisation of eigenvector
#--------------------------------------------------------------------
function normalization!(discretization::KSEDiscretization, 
                        U::AbstractArray{<:Real}, 
                        n::AbstractMatrix{<:Real})
    @unpack M₀ = discretization.matrices
    @unpack lₕ, nₕ = discretization
    @inbounds for k ∈ 1:nₕ
        @inbounds for σ ∈ 1:2
            @inbounds for l ∈ 1:lₕ+1   
                if !iszero(n[l,k,σ])
                    @views Ulkσ = U[:,k,l,σ]
                    coeff = sqrt(Ulkσ'*M₀*Ulkσ)
                    Ulkσ .= Ulkσ .* 1.0/coeff
                end
            end
        end
    end
    nothing
end


function normalization!(discretization::KSEDiscretization, 
                        U::AbstractArray{<:Real}, 
                        l::Int, 
                        k::Int,
                        σ::Int = 1)
    @unpack M₀ = discretization.matrices
    @views Ulkσ = U[:,k,l+1,σ]
    coeff = sqrt(Ulkσ'*M₀*Ulkσ)
    Ulkσ .= Ulkσ .* 1.0/coeff
    nothing
end


#####################################################################
#                          TYPES ORBITALS
#####################################################################

multiplicity(::KSEDiscretization, l::Int) = 4l+2

function orbitals_repartion(kd::KSEDiscretization, n::AbstractMatrix{<:Real})
    @unpack lₕ, nₕ, Nₕ = kd
    repart_orbitals = zeros(Int, 3, lₕ+1)
    for l ∈ 1:lₕ+1
        Nf = 0
        Np = 0
        mult = multiplicity(kd, l-1)
        for k ∈ 1:nₕ
            if n[l, k] == mult
                Nf += 1
            elseif n[l,k] > 0
                Np += 1
            end
        end
        repart_orbitals[1,i] = Nf
        repart_orbitals[2,i] = Np
        repart_orbitals[3,i] = Nₕ - Nf - Np
    end
    repart_orbitals
end