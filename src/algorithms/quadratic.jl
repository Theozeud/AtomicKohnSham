#--------------------------------------------------------------------
#                          Quadratic Method
#--------------------------------------------------------------------
"""
    QuadraticMethod

Implement the quadratic algorithm of "Quadratically convergent algorithm for fractional occupation numbers
in density functional theory" by Cancès, Kudin, Scuseria, Turinici. 
"""
struct QuadraticMethod <: SCFAlgorithm end

name(::QuadraticMethod) = "Quadratic-CancesTurinici"

"""
    CacheQuadratic

Cache structure to handle the `QuadraticMethod`.
"""
struct CacheQuadratic{densityType <: AbstractArray{<:Real},
    densitymatrixType <: AbstractArray{<:Real},
    operatorType <: AbstractArray{<:Real},
    singopType <: AbstractMatrix{<:Real},
    ΛType <: AbstractArray,
    T <: Real} <: SCFCache

    #= 
    DATA TO HANDLE THE UPDATE OF THE DENSITY MATRIX

    Let us recal that we have 
                    Γmo = (m * I    0  0
                             0      Λ  0
                             0      0  0 )   
    where we have :

            - Γmo   : the density matrix in the molecular basis
                                Γmo = Ω' * Γao * Ω
            - m     : the multiplicty of orbitals.

    The size of Γmo is (Nf + Np + Nv) x (Nf + Np + Nv)
    where :
            - Nf    : Number of fully occupied orbitals,
            - Np    : Number of partially occupied orbitals,
            - Nv    : Number of virtual orbitals.

    At each iteration we do : Ω = Ω * exp(A) and Λ = Λ + M
    To compute A and M, you have to solve the linear system L(A,M) = B
    in a vectorize way : B -> vecB which gives vecX -> (A,M)
    =#

    # DIAGONALIZED DENSITY MATRIX :  Γmo = Ω' * Γao * Ω 
    Γao::densitymatrixType
    Γmo::Any
    Ω::operatorType
    A::operatorType
    Λ::ΛType
    M::ΛType

    B::operatorType
    m::Int

    # REPARTITION ORBITALS
    repart_orbitals::Matrix{Int}    # Number of types of orbital for each section 
    # of the density matrix Γ :
    #   - Nf : number of fully occupied orbitals,
    #   - Np : number of partially occupied orbitals,
    #   - Nv : number of virtual orbitals. 
    # VECTOR CONVERSION FOR THE LINEAR SYSTEM
    vecA::Vector{T}
    vecB::Vector{T}
    indexconv::Vector{Tuple{Int, Int}}
end

function slice_orbitals(Nf::Int, Np::Int, Nv::Int)
    (1:Nf, (Nf + 1):(Nf + Np), (Nf + Np + 1):(Nf + Np + Nv))
end

function vec_slice_orbitals(Nf::Int, Np::Int, Nv::Int)
    Nvf = Nv * Nf
    Nvp = Nv * Np
    Npf = Np * Nf
    Npp = Np * Np
    (1:Nvf,
        (Nvf + 1):(Nvf + Nvp),
        (Nvf + Nvp + 1):(Nvf + Nvp + Npf),
        (Nvf + Nvp + Npf + 1):(Nvf + Nvp + Npf + Npp))
end

struct QuadraticOp
    Fmo::opType

    function QuadraticOp()
    end
end

function (Op::QuadraticOp)(X)
    Avf, Avp, Apf, M = X

    @unpack _commutator!(tmp1, A, Dmo) = Op
end

#--------------------------------------------------------------------
#                        Initialization
#--------------------------------------------------------------------

function create_cache_alg(::Quadratic,
        discretization::KohnShamDiscretization,
        rcacache::RCACache)
    @unpack matrices, lₕ, elT = discretization
    @unpack n, Noccup = rcacache
    @unpack M₀, S, H = matrices

    Nf, Np, Nv = Noccup

    # ORBITALS NUMBERS
    repart_orbitals = orbitals_repartion(discretization, n)

    # VECTOR 
    Avf = zeros(elT, Nv, Nf)
    Avp = zeros(elT, Nv, Np)
    Apf = zeros(elT, Np, Nf)
    M = zeros(elT, Np, Np)
    vecA = [Avf, Avp, Apf, M]

    vecB = zeros(elT, i)

    # DENSITY MATRIX
    Γao = zero_density_matrix(discretization)

    density_matrix!(discretization, U, n, Γ)

    m = multiplicty(discretization)
    Λ = Vector{Matrix{elT}}(undef, nbsection(discretization))

    for I in eachsection(discretization)
        @views ΓaoI = Γao[:, :, I...]
        @views ΩI = Ω[:, :, I...]
        @views ΛI = Λ[I]
        SΓS = Sq*ΓaoI*Sq
        Δ, V = eigen(SΓS)
        ΩI .= S * V[:, end:-1:1]
        Λ = diagm(Δ[(Nv + 1):(Nv + Np)])
        push!(Λ, ΛI)
    end
    M = zero(Λ)

    # LINEAR SYSTEM
    A = init_operator(discretization)
    B = init_operator(discretization)

    CacheQuadratic{typeof(D),
        typeof(Γao),
        typeof(Jd),
        typeof(Λ),
        elT}(
        D,
        Γao,
        Ω,
        A,
        Λ,
        M,
        B,
        m,
        repart_orbitals,
        vecX,
        vecB,
        indexconv)
end

#--------------------------------------------------------------------
#                             PERFORMSTEP
#--------------------------------------------------------------------

function performstep!(cache::CacheQuadratic, ::Quadratic, solver::KohnShamSolver)
    @unpack discretization, model, opts, energies = solver

    # STEP 1 : PREPARE THE LINEAR SYSTEM
    prepare_linear_system!(cache, discretization, model, opts.hartree)

    # STEP 2 : SOLVE THE LINEAR SYSTEM
    solve_linear_system!(cache, method)

    # STEP 3 : COMPUTE NEW Ω AND M
    update_density!(cache, discretization)

    # STEP 4 : UPDATE ENERGY
    update_energy!(solver)
end

function prepare_linear_system!(cache::CacheQuadratic,
        discretization::KohnShamDiscretization,
        model::KohnShamExtended,
        hartree::Real = 1.0)
    @unpack quadop, Λ, m, Γmo, Γao, Ω, Fmo, B = cache
    @unpack tmp1, tmp2 = quadop

    # COMPUTE Γmo 
    block(Γmo)[2] .= Λ

    # COMPUTE Γao
    _mul!(Γao, Ω, Γmo, Ω', tmp1)

    # COMPUTE Fao
    prepare_eigenvalue_problem!(discretization, model, D, hartree)
    @. Fao = discretization.matrices.H

    # COMPUTE Fmo
    _mul!(Fmo, Ω', H, Ω, tmp1)

    # COMPUTE B 
    _commutator!(tmp2, Fmo, Γmo, tmp1)

    @views Bpf = tmp2
    B[3]

    nothing
end

function linearsystem(cache::CacheQuadratic,
        discretization::KohnShamDiscretization,
        vecIN::AbstractVector)
    @unpack Γmo, Ω, A, M, Fmo = cache

    @unpack tmp1, tmp2, d, Jd, Q, tmpX, vecOUT = cache

    for I in eachsection(discretization)
        @views Av = A[:, :, I...]
        @views Ωv = Ω[:, :, I...]
        @views dv = d[:, :, I...]
        @views Γv = Γmo[:, :, I...]
        @views Mv = M[I]

        Nf, Np, Nv = Noccup[:, I...]
        _, slicep, _ = slice_orbitals(Nf, Np, Nv)

        # Take vecINv

        # CONVERT X INTO (A,M)
        vecIN_to_AM!(Av, Mv, vecINv, Nf, Np, Nv)

        # COMPUTE Q(M)
        remove_trace!(Mv)

        # COMPUTE d = Ω × ([A,ΓM] + Q(M)) × Ω'
        _commutator!(tmp2, Av, Γv, tmp1)
        @views tmp2pp = tmp2[slicep, slicep]
        @. tmp2pp += Mv
        _mul!(dv, Ωv, tmp2, Ωv', tmp1)
    end

    # COMPUTE THE TRIAL DENSITY d
    density!(discretization, d, rhod)

    # COMPUTE HARTREE AND EXCHANGE CORRELATION FOR D
    # Jd, Qd

    for I in eachsection(discretization)
        @views Av = A[:, :, I...]
        @views Ωv = Ω[:, :, I...]
        @views Fmov = Fmo[:, :, I...]
        @views Γv = Γmo[:, :, I...]
        @views Mv = M[I]

        Nf, Np, Nv = Noccup[:, I...]
        _, slicep, _ = slice_orbitals(Nf, Np, Nv)

        # COMPUTE Gmo = Ω' × (J(d) + Q_xc(DM)) × Ω
        tmp1 .= Jd .+ Qd
        _mul!(Gmo, Ωv', tmp1, Ωv, tmp2)

        # COMPUTE Z = [Fmo,A] + GM
        _commutator!(tmpZ, Fmo, Av, tmp1)
        mul!(tmp1, tmpZ, Γmo)
        @. tmpZ += GM

        # COMPUTE   X = 0.5 × ([Fmo, A]Γmo + Fmo[A,Γmo]) + GM×ΓM + FM×Q(M) 
        _commutator!(tmp2, A, ΓM, tmp1)
        mul!(tmpX, FM, tmp2)
        @. tmpX += tmp1
        @. tmpX *= 1/2
        mul!(tmp1, GM, ΓM)
        @. tmpX += tmp1
        @views tmp1_p = tmp1[:, slicep]
        @views FM_p = FM[:, slicep]
        @views tmpX_p = tmpX[:, slicep]
        mul!(tmp1_v, FM_v, Mv)
        @. tmpX_v += tmp1_v

        # OUTPUT
        tmp1 .= tmpX .- tmpX'
        @views Zpp = Z[slicep, slicep]
        remove_trace!(Zpp)
        AM_to_vecOUT!(vecOUT, tmp1, Zpp, Nf, Np, Nv)
    end

    vecOUT
end

function solve_linear_system!(cache::CacheQuadratic, method::Quadratic)
    @unpack vecIn, vecB, A, M, Nf, Np, Nv = cache

    sol = linsolve(quadop, B; issymetric = true, isposdef = true)

    vecIn .= sol[1]

    for I in eachsection(discretization)
        @views Av = A[:, :, I...]
        @views Mv = M[I]
        Nf, Np, Nv = Noccup[:, I...]
        vecIn_to_AM!(A, M, vecX, Nf, Np, Nv)
    end

    nothing
end

function update_density!(cache::CacheQuadratic)
    @unpack Ω, Λ, M, A = cache

    # UPDATE Ω
    Ω .= Ω * exp(A)
    Λ
    # UPDATE Λ 
    block(Γv)[2] .+= M

    @unpack quadop, Λ, m, Γmo, Γao, Ω, Fmo, B = cache
    @unpack tmp1, tmp2 = quadop

    # COMPUTE Γmo 
    block(Γmo)[2] .= Λ

    # COMPUTE Γao
    _mul!(Γao, Ω, Γmo, Ω', tmp1)
end

#####################################################################
#          CONVERSION BETWEEN MATRIX AND VECTOR REPRESENTATIONS
#####################################################################

function X_to_AM!(
        A::AbstractMatrix, M::AbstractMatrix, X::AbstractVector, Nf::Int, Np::Int, Nv::Int)
    slicef, slicep, slicev = slice_orbitals(Nf, Np, Nv)
    slicevf, slicevp, slicepf, slicepp = vec_slice_orbitals(Nf, Np, Nv)

    _copy_vec_to_mat!(A, slicev, slicef, X, slicevf)
    _copy_vec_to_mat!(A, slicev, slicep, X, slicevp)
    _copy_vec_to_mat!(A, slicep, slicef, X, slicepf)
    _copy_vec_to_mat!(M, axes(M, 1), axes(M, 2), X, slicepp)

    LinearAlgebra._copy_adjtrans!(A, slicef, slicev, A, slicev, slicef, antiadjoint)
    LinearAlgebra._copy_adjtrans!(A, slicep, slicev, A, slicev, slicep, antiadjoint)
    LinearAlgebra._copy_adjtrans!(A, slicef, slicep, A, slicep, slicef, antiadjoint)
    nothing
end

function AM_to_X!(
        X::AbstractVector, A::AbstractMatrix, M::AbstractMatrix, Nf::Int, Np::Int, Nv::Int)
    slicef, slicep, slicev = slice_orbitals(Nf, Np, Nv)
    slicevf, slicevp, slicepf, slicepp = vec_slice_orbitals(Nf, Np, Nv)

    _copy_mat_to_vec!(X, slicevf, A, slicev, slicef)
    _copy_mat_to_vec!(X, slicevp, A, slicev, slicep)
    _copy_mat_to_vec!(X, slicepf, A, slicep, slicef)
    _copy_mat_to_vec!(X, slicepp, M, axes(M, 1), axes(M, 2))
    nothing
end
