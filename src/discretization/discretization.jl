#--------------------------------------------------------------------
#                          KSE MATRICES
#--------------------------------------------------------------------
struct KSEMatrices{T <: Real,
    typeMatrix <: AbstractMatrix{<:Real},
    typeVectorMatrix <: Vector{<:AbstractMatrix{<:Real}},
    typeHam <: AbstractArray{<:Real}}
    # FEM MATRICES
    A::typeMatrix                   # Matrix of Qᵢ'Qⱼ'
    M₀::Matrix{T}                   # Matrix of QᵢQⱼ
    M₋₁::typeMatrix                 # Matrix of 1/x QᵢQⱼ
    M₋₂::typeMatrix                 # Matrix of 1/x² QᵢQⱼ
    F::Dict{Tuple{Int, Int, Int}, T}   # Tensor of 1/x QᵢQⱼQₖ
    S::Matrix{T}                    # √(M₀⁻¹)
    Sinv::Matrix{T}                 # √(M₀)
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
struct DiscretizationCache{T <: Real, Tρ}
    tmp_MV::Matrix{T}               # Store the contraction F:C
    tmp_B::Vector{T}                # Matrix of 4πρQᵢ
    tmp_C::Vector{T}                # Matrix for V(ρ)
    tmp_vect::Vector{T}             # Vector to store temporary vector
    tmp_ρ::Tρ                        # Matrix of one or two rows (depending on the spin
    # polarization to store density evaluations during
    # quadrature method.
    #   First row  -> ρ↑
    #   Second row -> ρ↓
    tmp_vρ::Tρ                      # Similary to tmpρ but to store the exchange-correlation
    # potential.
    tmp_vρ2::Tρ
end

function create_cache_discretization(lₕ::Int, Nₕ::Int, T::Type, n_spin::Int, n_eval::Int)
    # FEM MATRICES
    A = spzeros(T, Nₕ, Nₕ)
    M₀ = zeros(T, Nₕ, Nₕ)
    M₋₁ = spzeros(T, Nₕ, Nₕ)
    M₋₂ = spzeros(T, Nₕ, Nₕ)
    F = Dict{Tuple{Int, Int, Int}, T}()
    S = zeros(T, Nₕ, Nₕ)
    Sinv = zeros(T, Nₕ, Nₕ)
    # MATRICES COMPOSING THE HAMILTONIAN
    H = flexible_zeros(T, (Nₕ, Nₕ, lₕ+1), n_spin)
    Kin = _spzeros(T, Nₕ, Nₕ, lₕ+1)
    Coulomb = spzeros(T, Nₕ, Nₕ)
    Hfix = _spzeros(T, Nₕ, Nₕ, lₕ+1)
    Hartree = spzeros(T, Nₕ, Nₕ)
    VxcUP = spzeros(T, Nₕ, Nₕ)
    VxcDOWN = spzeros(T, Nₕ, Nₕ)
    # Initialization of array for temporary stockage of computations
    tmp_MV = zeros(T, Nₕ, Nₕ)
    tmp_B = zeros(T, Nₕ)
    tmp_C = zeros(T, Nₕ)
    tmp_vect = zeros(T, Nₕ)
    tmp_ρ = flexible_zeros(T, n_spin, (n_eval,))
    tmp_vρ = flexible_zeros(T, n_spin, (n_eval,))
    tmp_vρ2 = flexible_zeros(T, n_spin, (n_eval,))
    KSEMatrices{T, typeof(A), typeof(Hfix), typeof(H)}(A, M₀, M₋₁, M₋₂, F, S, Sinv,
        H, Kin, Coulomb, Hfix, Hartree, VxcUP, VxcDOWN),
    DiscretizationCache{T, typeof(tmp_ρ)}(
        tmp_MV, tmp_B, tmp_C, tmp_vect, tmp_ρ, tmp_vρ, tmp_vρ2)
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
    M <: Mesh,
    Mat <: KSEMatrices,
    FIM <: FEMIntegrationMethod}
    lₕ::Int
    nₕ::Int
    Nₕ::Int
    basis::B
    mesh::M
    Rmin::T
    Rmax::T
    n_spin::Int
    matrices::Mat
    cache::DiscretizationCache{T}
    fem_integration_method::FIM     # Integration method to compute integrals
    # for fem's matrices
    function KSEDiscretization(lₕ::Int,
            basis::FEMBasis,
            mesh::Mesh,
            n_spin::Int,
            nₕ::Int = length(basis),
            fem_integration_method::FEMIntegrationMethod = GaussLegendre(basis))
        elT = eltype(basis)
        Nₕ = length(basis)
        n_eval = fem_integration_method.npoints
        matrices, cache = create_cache_discretization(lₕ, Nₕ, elT, n_spin, n_eval)
        new{elT,
            typeof(basis),
            typeof(mesh),
            typeof(matrices),
            typeof(fem_integration_method)}(lₕ,
            nₕ,
            Nₕ,
            basis,
            mesh,
            elT(first(mesh)),
            elT(last(mesh)),
            n_spin,
            matrices,
            cache,
            fem_integration_method)
    end
end

function init_cache!(discretization::KSEDiscretization, model::KSEModel)
    @unpack lₕ, basis, matrices = discretization
    @unpack A, M₀, M₋₁, M₋₂, F, S, Sinv, Kin, Coulomb, Hfix = matrices

    # CREATION OF FEM MATRICES
    fill_stiffness_matrix!(basis, A)
    fill_mass_matrix!(basis, M₀)
    fill_mass_matrix!(basis, -1, M₋₁)
    lₕ == 0 || fill_mass_matrix!(basis, -2, M₋₂)
    iszero(model.hartree) || fill_mass_tensor!(basis, -1, F)
    S .= sqrt(inv(Symmetric(M₀)))
    Sinv .= sqrt(Symmetric(M₀))

    # CREATION OF THE FIX PART OF THE HAMILTONIAN
    kinetic_matrix!(discretization)
    coulomb_matrix!(discretization, model)
    for l in 1:(lₕ + 1)
        @. Hfix[l] = Kin[l] + Coulomb
    end
    nothing
end

#--------------------------------------------------------------------
#                               API
#--------------------------------------------------------------------

Base.eltype(::KSEDiscretization{T}) where {T} = T
dim(discretization::KSEDiscretization) = discretization.Nₕ * (discretization.lₕ + 1)
multiplicity(::KSEDiscretization) = 2

zero_density(kd::KSEDiscretization) = flexible_zeros(eltype(kd), (kd.Nₕ, kd.Nₕ), kd.n_spin)
function zero_orbitals(kd::KSEDiscretization)
    flexible_zeros(eltype(kd), (kd.Nₕ, kd.nₕ, kd.lₕ+1), kd.n_spin)
end
function zero_orbitals_energy(kd::KSEDiscretization)
    flexible_zeros(eltype(kd), (kd.lₕ+1, kd.nₕ), kd.n_spin)
end
function zero_occupation_number(kd::KSEDiscretization)
    flexible_zeros(eltype(kd), (kd.lₕ+1, kd.nₕ), kd.n_spin)
end
function zero_density_matrix(kd::KSEDiscretization)
    flexible_zeros(eltype(kd), (kd.Nₕ, kd.Nₕ, kd.lₕ+1), kd.n_spin)
end

function zero_energies(kd::KSEDiscretization, model::KSEModel)
    elT = eltype(kd)
    d = Dict(:Etot => zero(elT),                                     # Total energy
        :Ekin => zero(elT),                                     # Kinetic energy
        :Ecou => zero(elT),                                     # Coulomb energy
        :Ehar => zero(elT))                                     # Hartree energy
    if has_exchcorr(model)
        (d[:Eexc] = zero(elT))                                          # Exchange-correlation energy
        kd.n_spin != :lsda || (d[:Ekincorr] = zero(elT))                # Kinetic-correlation energy
    end
    d
end

zero_operator(kd::KSEDiscretization) = zeros(eltype(kd), kd.Nₕ, kd.Nₕ, kd.lₕ+1)
zero_single_operator(kd::KSEDiscretization) = zeros(eltype(kd), kd.Nₕ, kd.Nₕ)

#--------------------------------------------------------------------
#               Find Orbital : Solve the eigen problems
#--------------------------------------------------------------------
function build_hamiltonian!(discretization::KSEDiscretization,
        model::KSEModel,
        D::AbstractArray{<:Real},
        hartree::Real = true)
    @unpack H, Hfix, Hartree, VxcUP, VxcDOWN = discretization.matrices

    # COMPUTE HARTREE MATRIX
    iszero(hartree) || hartree_matrix!(discretization, D, hartree)

    # COMPUTE EXCHANGE CORRELATION MATRIX
    !has_exchcorr(model) || exchange_corr_matrix!(discretization, model, D)

    # BUILD THE HAMILTONIAN OF THE lᵗʰ SECTION FOR EACH SPIN σ ∈ {↑,↓}
    @threads for l in 0:discretization.lₕ
        @views vH = H[:, :, l + 1, 1]
        @. vH = Hfix[l + 1] + VxcUP + Hartree
    end

    if discretization.n_spin == 2
        @threads for l in 0:discretization.lₕ
            @views vH = H[:, :, l + 1, 2]
            @. vH = Hfix[l + 1] + VxcDOWN + Hartree
        end
    end
    nothing
end

function find_orbital!(discretization::KSEDiscretization,
        U::AbstractArray{<:Real},
        ϵ::AbstractArray{<:Real})
    @unpack lₕ, nₕ, matrices, n_spin = discretization
    @unpack S, H = matrices

    # SOLVE THE GENERALIZED EIGENVALUE PROBLEM FOR EACH SECTION l AND FOR EACH σ ∈ {↑,↓}
    for σ in 1:n_spin
        for l in 0:lₕ
            @views vH = H[:, :, l + 1, σ]
            λ, V = eigen(Symmetric(S*vH*S))
            @views ϵ[l + 1, :, σ] = λ[1:nₕ]
            @views U[:, :, l + 1, σ] = S*V[:, 1:nₕ]
        end
    end
end

#--------------------------------------------------------------------
#               Normamisation of eigenvector
#--------------------------------------------------------------------
function normalization!(discretization::KSEDiscretization,
        U::AbstractArray{<:Real},
        l::Int,
        k::Int,
        σ::Int = 1)
    @unpack M₀ = discretization.matrices
    @views Ulkσ = U[:, k, l + 1, σ]
    coeff = sqrt(Ulkσ'*M₀*Ulkσ)
    Ulkσ .= Ulkσ .* 1.0/coeff
    nothing
end

#####################################################################
#                          TYPES ORBITALS
#####################################################################

multiplicity(::KSEDiscretization, l::Int) = 4l+2

function orbitals_repartion(kd::KSEDiscretization, n::AbstractArray{<:Real})
    @unpack lₕ, nₕ, Nₕ = kd
    repart_orbitals = zeros(Int, 3, lₕ+1)
    for l in 1:(lₕ + 1)
        Nf = 0
        Np = 0
        mult = multiplicity(kd, l-1)
        for k in 1:nₕ
            if n[l, k] == mult
                Nf += 1
            elseif n[l, k] > 0
                Np += 1
            end
        end
        repart_orbitals[1, i] = Nf
        repart_orbitals[2, i] = Np
        repart_orbitals[3, i] = Nₕ - Nf - Np
    end
    repart_orbitals
end
