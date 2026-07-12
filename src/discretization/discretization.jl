# ===================================================================
#                          FEM OPERATORS
# ===================================================================
"""
    FEMOperators{T}

Finite element operators used to assemble the radial Kohn–Sham Hamiltonian.

This structure stores the discretized FEM operators associated with the radial
basis functions, including stiffness and mass matrices, as well as auxiliary
operators required for Coulomb and Hartree terms.

The operators are constructed once for a given basis and reused throughout the
self-consistent field (SCF) iterations.

# Fields
- `A`: Stiffness matrix ⟨∂r Qᵢ, ∂r Qⱼ⟩.
- `M₀`: Mass matrix ⟨Qᵢ, Qⱼ⟩.
- `M₋₁`: Weighted mass matrix ⟨Qᵢ, Qⱼ / r⟩.
- `M₋₂`: Weighted mass matrix ⟨Qᵢ, Qⱼ / r²⟩.
- `F`: Sparse tensor associated with ⟨Qᵢ Qⱼ Qₖ / r⟩.
- `S`: Square root of `M₀⁻¹`, used for orthonormalization.
- `Sinv`: Square root of `M₀`.

This structure is backend-agnostic and does not depend on the SCF algorithm.
"""
struct FEMOperators{T<:Real}
    A::SparseMatrixCSC{T,Int}
    M₀::Matrix{T}
    M₋₁::SparseMatrixCSC{T,Int}
    M₋₂::SparseMatrixCSC{T,Int}
    F::Dict{Tuple{Int, Int, Int}, T}
    S::Matrix{T}
    Sinv::Matrix{T}

    function FEMOperators(Nₕ::Int, ::Type{T}) where T
        A       = spzeros(T, Nₕ, Nₕ)
        M₀      = zeros(T, Nₕ, Nₕ)
        M₋₁     = spzeros(T, Nₕ, Nₕ)
        M₋₂     = spzeros(T, Nₕ, Nₕ)
        F       = Dict{Tuple{Int, Int, Int}, T}()
        S       = zeros(T, Nₕ, Nₕ)
        Sinv    = zeros(T, Nₕ, Nₕ)
        new{T}(A, M₀, M₋₁, M₋₂, F, S, Sinv)
    end
end

# ===================================================================
#                      HAMILTONIAN COMPONENTS
# ===================================================================
"""
    KSHamiltonian{T, HamType}

Discrete Kohn–Sham Hamiltonian and its components.

This structure stores the block-structured Kohn–Sham Hamiltonian for all angular
momentum channels and spin components, together with its physical contributions
(kinetic, Coulomb, Hartree, and exchange–correlation terms).

The Hamiltonian blocks are updated during the SCF procedure, while some components
(e.g. kinetic and nuclear Coulomb terms) are fixed after initialization.

# Fields
- `H`: Full Kohn–Sham Hamiltonian blocks indexed by angular momentum and spin.
- `Kin`: Kinetic energy operators for each angular momentum channel.
- `Coulomb`: Nuclear Coulomb potential operator.
- `Hfix`: Fixed part of the Hamiltonian (`Kin + Coulomb`).
- `Hartree`: Hartree potential operator.
- `VxcUP`: Exchange–correlation potential for spin-up electrons.
- `VxcDOWN`: Exchange–correlation potential for spin-down electrons.

When spin polarization is disabled, only `VxcUP` is used.
"""
struct KSHamiltonian{T<:Real, HamType}
    H::HamType
    Kin::Vector{SparseMatrixCSC{T, Int}}
    Coulomb::SparseMatrixCSC{T,Int}
    Hfix::Vector{SparseMatrixCSC{T, Int}}
    Hartree::SparseMatrixCSC{T,Int}
    VxcUP::SparseMatrixCSC{T,Int}
    VxcDOWN::SparseMatrixCSC{T,Int}

    function KSHamiltonian(Nₕ::Int, lₕ::Int, nspin::Int, ::Type{T}) where T
        H       = flexible_zeros(T, (Nₕ, Nₕ, lₕ+1), nspin)
        Kin     = _spzeros(T, Nₕ, Nₕ, lₕ+1)
        Coulomb = spzeros(T, Nₕ, Nₕ)
        Hfix    = _spzeros(T, Nₕ, Nₕ, lₕ+1)
        Hartree = spzeros(T, Nₕ, Nₕ)
        VxcUP   = spzeros(T, Nₕ, Nₕ)
        VxcDOWN = spzeros(T, Nₕ, Nₕ)
        new{T,typeof(H)}(H,Kin,Coulomb,Hfix, Hartree, VxcUP, VxcDOWN)
    end
end

# ===================================================================
#                           KSE CACHE
# ===================================================================
"""
Workspace for Hartree potential computations.
"""
struct HartreeWorkspace{T <: Real}
    W::Vector{T}                # Matrix for V(ρ)
    FW::Matrix{T}               # Store the contraction F:W
    B::Vector{T}                # Matrix of 4πρQᵢ

    function HartreeWorkspace(Nₕ::Int, ::Type{T}, hartree::Real) where T
        if iszero(hartree)
            W  = zeros(T, 0)
            FW = zeros(T, 0, 0)
            B  = zeros(T, 0)
            return new{T}(W,FW,B)
        else
            W  = zeros(T, Nₕ)
            FW = zeros(T, Nₕ, Nₕ)
            B  = zeros(T, Nₕ)
            return new{T}(W,FW,B)
        end
    end
end

"""
Workspace for exchange–correlation evaluations.
"""
struct ExcWorkspace{Tρ}
    # Matrix of one or two rows (depending on the spin
    # polarization to store density evaluations during
    # quadrature method.
    #   First row  -> ρ↑
    #   Second row -> ρ↓
    ρ_buf::Tρ
    # Similary to tmpρ but to store the exchange-correlation
    # potential.
    vρ_buf::Tρ
    vρ_buf2::Tρ

    function ExcWorkspace(neval::Int, nspin::Int, ::Type{T}, exc::Bool) where T
        if exc
            ρ_buf = flexible_zeros(T, nspin, (neval,))
            vρ_buf = flexible_zeros(T, nspin, (neval,))
            vρ_buf2 = flexible_zeros(T, nspin, (neval,))
            return new{typeof(ρ_buf)}(ρ_buf, vρ_buf, vρ_buf2)
        else
            ρ_buf = flexible_zeros(T, 0, (0,))
            vρ_buf = flexible_zeros(T, 0, (0,))
            vρ_buf2 = flexible_zeros(T, 0, (0,))
            return new{typeof(ρ_buf)}(ρ_buf, vρ_buf, vρ_buf2)
        end
    end
end

"""
Temporary buffers for numerical evaluations.
"""
struct EvalWorkSpace{T}
    buf1::Vector{T}
    buf2::Vector{T}

    function EvalWorkSpace(Nₕ::Int, ::Type{T}) where T
        buf1  = zeros(T, Nₕ)
        buf2  = zeros(T, Nₕ)
        new{T}(buf1, buf2)
    end
end

"""
    DiscretizationCache

Internal workspaces for discretization-dependent computations.

This structure groups all temporary buffers required to assemble the Hartree and
exchange–correlation contributions and to evaluate densities and potentials during
SCF iterations.

The cache is fully preallocated at construction time in order to avoid memory
allocations inside performance-critical loops.

This structure is internal and should not be modified by users.
"""
struct DiscretizationCache{T <: Real, Tρ}
    hartw::HartreeWorkspace{T}
    excw::ExcWorkspace{Tρ}
    evalw::EvalWorkSpace{T}

    function DiscretizationCache(Nₕ::Int, T::Type, nspin::Int, neval::Int, hartree::Real,
                                exc::Bool)
        # Worspace for Hartree computations
        hartw   = HartreeWorkspace(Nₕ, T, hartree)
        # Workspace for exchange-correlation computations
        excw    = ExcWorkspace(neval, nspin,T, exc)
        # Workspace for evaluationS
        evalw   = EvalWorkSpace(Nₕ, T)
        new{T, typeof(excw.ρ_buf)}(hartw, excw, evalw)
    end
end

# ===================================================================
#                          KSE Discretization
# ===================================================================
"""
    KSEDiscretization

Finite element discretization of the radial Kohn–Sham equations with spherical
symmetry.

This structure defines the numerical discretization used to represent radial
Kohn–Sham orbitals in a basis of finite element functions combined with spherical
harmonics. It stores the FEM operators, Kohn–Sham Hamiltonian components, and
associated workspaces required for SCF computations.

The discretization is independent of the SCF algorithm and can be reused with
different solvers.

# Mathematical form

The radial part of a Kohn–Sham orbital is expanded as

    u(r,θ,φ) = ∑_{l,m,n} u_{nlm} Qₙ(r) / r · Yₗᵐ(θ,φ),

where `(Qₙ)` are radial FEM basis functions and `Yₗᵐ` are spherical harmonics.

# Fields
- `lₕ`: Maximum angular momentum quantum number.
- `nₕ`: Number of orbitals per angular momentum channel.
- `Nₕ`: Number of radial basis functions.
- `nspin`: Number of spin components.
- `N`: Total number of electrons.
- `Rmax`: Radial cutoff.
- `basis`: Radial FEM basis.
- `femops`: Finite element operators.
- `ksham`: Kohn–Sham Hamiltonian components.
- `cache`: Preallocated workspaces.
- `fem_integration_method`: Quadrature method for FEM integrals.
- `eigensolver`: Strategy used to solve the Kohn–Sham eigenproblem in each
  `(l, σ)` channel — [`PartialEigenSolver`](@ref) (default) or
  [`FullEigenSolver`](@ref).

Use `init_cache!` to assemble all discretization-dependent operators.
"""
struct KSEDiscretization{T <: Real, B <: FEMBasis, O <: FEMOperators, H <: KSHamiltonian,
                        FIM <: FEMIntegrationMethod, ES <: EigenSolver,}

    lₕ::Int
    nₕ::Int
    Nₕ::Int
    nspin::Int
    N::T
    Rmax::T
    basis::B
    femops::O
    ksham::H
    cache::DiscretizationCache{T}
    fem_integration_method::FIM
    eigensolver::ES
    function KSEDiscretization(basis::FEMBasis, model::KSEModel; lh::Int, nh::Int = 10,
                    fem_integration_method::FEMIntegrationMethod = GaussLegendre(basis),
                    eigensolver::EigenSolver = PartialEigenSolver())
        @unpack nspin, hartree, N = model
        T = eltype(basis)
        Nₕ = length(basis)
        Rmax = T(last(basis.mesh))

        femops = FEMOperators(Nₕ, T)
        ksham = KSHamiltonian(Nₕ, lh, nspin, T)

        neval = fem_integration_method.npoints
        exc = has_exchcorr(model)
        cache = DiscretizationCache(Nₕ, T, nspin, neval, hartree, exc)

        new{T,
            typeof(basis),
            typeof(femops),
            typeof(ksham),
            typeof(fem_integration_method),
            typeof(eigensolver)}(lh, nh, Nₕ, nspin, N, Rmax, basis, femops,
            ksham, cache, fem_integration_method, eigensolver)
    end
end


"""
    init_cache!(discretization, model)

Assemble all finite element operators and initialize the fixed components of the
Kohn–Sham Hamiltonian for the given model.

This function must be called once before starting the SCF iterations.
"""
function init_cache!(discretization::KSEDiscretization, model::KSEModel)
    @unpack lₕ, basis, femops, ksham = discretization
    @unpack A, M₀, M₋₁, M₋₂, F, S, Sinv = femops
    @unpack Kin, Coulomb, Hfix = ksham

    # Creation of the fem operators
    fill_stiffness_matrix!(basis, A)
    fill_mass_matrix!(basis, M₀)
    fill_mass_matrix!(basis, -1, M₋₁)
    lₕ == 0 || fill_mass_matrix!(basis, -2, M₋₂)
    iszero(model.hartree) || fill_mass_tensor!(basis, -1, F)
    S .= sqrt(inv(Symmetric(M₀)))
    Sinv .= sqrt(Symmetric(M₀))

    # Creation of the fix part of the hamiltonian
    assemble_kinetic!(discretization)
    assemble_coulomb!(discretization, model.Z)
    for l in 1:(lₕ + 1)
        @. Hfix[l] = Kin[l] + Coulomb
    end
    nothing
end

# ===================================================================
#                               API
# ===================================================================
Base.eltype(::KSEDiscretization{T}) where {T} = T

"""
Allocate a zero-initialized density matrix compatible with the discretization.
"""
function zero_density(kd::KSEDiscretization)
    flexible_zeros(eltype(kd), (kd.Nₕ, kd.Nₕ), kd.nspin)
end

"""
Allocate zero-initialized Kohn–Sham orbital coefficients.
"""
function zero_orbitals_coeffs(kd::KSEDiscretization)
    flexible_zeros(eltype(kd), (kd.Nₕ, kd.nₕ, kd.lₕ+1), kd.nspin)
end

"""
Allocate zero-initialized Kohn–Sham orbital energies.
"""
function zero_orbitals_energies(kd::KSEDiscretization)
    flexible_zeros(eltype(kd), (kd.lₕ+1, kd.nₕ), kd.nspin)
end

"""
Allocate zero-initialized occupation numbers for Kohn–Sham orbitals.
"""
function zero_occupation_numbers(kd::KSEDiscretization)
    flexible_zeros(eltype(kd), (kd.lₕ+1, kd.nₕ), kd.nspin)
end

"""
Convert a flattened orbital index into angular momentum, radial, and spin indices.
"""
function convert_index(discretization::KSEDiscretization, idx::Int)
    @unpack lₕ, nₕ, nspin = discretization
    if nspin == 1
        l = rem(idx - 1, lₕ+1)
        k = div(idx-1, lₕ+1)+1
        return (l, k, 1)
    else
        σ = div(idx-1, (lₕ+1)*nₕ)+1
        idxσ = rem(idx-1, (lₕ+1)*nₕ)
        l = rem(idxσ, lₕ+1)
        k = div(idxσ, lₕ+1)+1
        return (l, k, σ)
    end
end

"""
    convert_index_nl(discretization, idx)

Map the flattened orbital index `idx` to quantum numbers.

Returns `(n, l)` for spin-unpolarized discretizations and `(n, l, σ)` for
spin-polarized ones, with the convention `n > l` and `n = k + l` (where `k ≥ 1`
is the radial index).
"""
function convert_index_nl(discretization::KSEDiscretization, idx::Int)
    @unpack lₕ, nₕ, nspin = discretization

    if nspin == 1
        l = rem(idx - 1, lₕ + 1)
        k = div(idx - 1, lₕ + 1) + 1
        n = k + l
        return (n, l)
    else
        σ    = div(idx - 1, (lₕ + 1) * nₕ) + 1
        idxσ = rem(idx - 1, (lₕ + 1) * nₕ)
        l    = rem(idxσ, lₕ + 1)
        k    = div(idxσ, lₕ + 1) + 1
        n    = k + l
        return (n, l, σ)
    end
end

"""
Return the (2l+1) or (4l+2) degeneracy associated with a given orbital.
"""
function degeneracy(discretization::KSEDiscretization, idx::Int)
    @unpack nspin = discretization
    l, _ = convert_index(discretization, idx)
    if nspin == 1
        return 4 * l + 2
    else
        return 2 * l + 1
    end
end


#=
# ===================================================================
#               Normamisation of eigenvector
# ===================================================================
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
=#
