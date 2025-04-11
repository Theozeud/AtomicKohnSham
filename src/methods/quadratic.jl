#####################################################################
#                          Quadratic Method
#####################################################################

struct Quadratic{typeMethod, typeOpts<:NamedTuple} <: SCFMethod 
    method::typeMethod  # Method  to solve the linear system
    opts::typeOpts      # Options to solve the linear system
end

name(::Quadratic) = "Quadratic" 


struct CacheQuadratic{densityType <: AbstractArray} <: SCFCache

    Nf::Int             # Number of fully occupied orbitals
    Np::Int             # Number of partially occupied orbitals
    Nv::Int             # Number of virtual orbitals 

    # DIAGONALIZED DENSITY MATRIX :  Γdiag = Ω' * Γ * Ω
    Γdiag::densityType     
    Ω::densityType       
    
    # DATA TO HANDLE THE UPDATE OF THE DENSITY MATRIX
    # Let us recal that we have 
    #                       Γdiag = (m * I    0  0
    #                                  0      Λ  0
    #                                  0      0  0 )
    # where m is the multiplicty of orbitals depending on the discretization.
    # At each iteration we do : Ω = Ω * exp(A) and Λ = Λ + M
    # To compute A and M, you have to solve the linear system L(A,M) = B.
    A::densityType
    M::densityType         
    B::densityType
    vecX::vectorType
    vecB::vectorType    
    
    FA::densityType
    FM::densityType
    
    # TEMPORARY MATRICES TO STORE INTERMEDIATE CALCULATIONS
    tmp1::densityType
    tmp2::densityType 
    tmpd::densityType  
    tmpX::densityType
    tmpQM::densityType
end


slice_orbitals(Nf::Int, Np::Int, Nv::Int) = (1:Nf, Nf+1:Nf+Np, Nf+Np+1:Nf+Np+Nv)

function vec_slice_orbitals(Nf::Int, Np::Int, Nv::Int)
    Nvf = Nv * Nf
    Nvp = Nv * Np
    Npf = Np * Nf
    Npp = Np * Np
    (1:Nvf, 
     Nvf+1:Nvf+Nvp, 
     Nvf+Nvp+1:Nvf+Nvp+Npf, 
     Nvf+Nvp+Npf+1:Nvf+Nvp+Npf+Npp)
end

#####################################################################
#                        Initialization
#####################################################################

function create_cache_method(method::Quadratic, discretization::KohnShamDiscretization)

    @unpack M₀, elT = discreztization

    # INIT OCCUPATIONS COUNT
    Nf, Np, Nv = Noccup

    # COMPUTE Γ
    Γ = init_density_matrix(discretization)
    density_matrix!(discretization, U, n, Γ)

    # INITIALIZATION OF Ω AND Λ 

    # 1. Compute the square root of the overlap matrix M₀
    S = sqrt(M₀)

    # 2. Construct SΓS by applying S * Γ * S block by block
    SΓS = similar(Γ)
    for i in 1:nblocks(Γ)
        @views Γi = block(Γ)[i]
        @views SΓSi = block(SΓS)[i]
        SΓSi .= S * Γi * S
    end

    # 3. Perform eigen decomposition of SΓS
    Δ, V = eigen(SΓS)

    # 4. Construct Ω and Λ
    Ω = inv(S) * V[:, end:-1:1]  
    Λ = diagm(Δ[Nv+1:Nv+Np])     

    # 5. Construct Γ_decomp as a block diagonal matrix
    Γdiag = BlockDiagonal([
                2 * Eye(Nf),         # Block of size (Nf, Nf)
                Λ,                   # Diagonal block (Np, Np)
                Zeros(elT, Nv, Nv)   # Zero block (Nv, Nv)
            ])

    # INIT "CACHES VARIABLES
    A = zero(D)
    M = zero(Λ)
    B = 

    CacheQuadratic{typeof(D)}(  Nf, Np, Nv, Γdiag, Ω, A, M, B, vecX, vecB, FAO, 
                                tmp1, tmp2, tmpd, tmpX, tmpQM)
end

#####################################################################
#                        Performstep
#####################################################################

function performstep!(cache::CacheQuadratic, ::Quadratic, solver::KohnShamSolver)
    # STEP 1 : PREPARE THE LINEAR SYSTEM
    prepare_linear_system!(cache)

    # STEP 2 : SOLVE THE LINEAR SYSTEM
    solve_linear_system!(cache, method)

    # STEP 3 : COMPUTE NEW Ω AND M
    update_density!(cache)

    # STEP 4 : UPDATE ENERGY
    update_energy!(solver)
end

#####################################################################
#                        LINEAR SYSTEM
#####################################################################
function prepare_linear_system!(cache::CacheQuadratic)
    @unpack Nf, Np, Nv,
            ΓM, Ω, A, M, B, FA, FM  = cache

    @unpack tmp1 = cache

    _, slicep, _ = slice_orbitals(Nf, Np, Nv)

    # COMPUTE FA


    # COMPUTE FM
    _mul!(FM, Ω', FA, Ω, tmp1)

    # COMPUTE B
    _commutator!(B, FM, ΓM, tmp1)
    @views Bpp = B[slicep, slicep]
    remove_trace!(Bpp)
    _copy_mat_to_vec!(vecB, axes(vecB,1), B, axes(B,1), axes(B,2))

end


function linearsystem(  cache::CacheQuadratic, 
                        X::AbstractVector)

    @unpack Nf, Np, Nv,
            ΓM, Ω, A, M, FM  = cache

    @unpack tmp1, tmp2, tmpQM, tmpd, tmpX, vecX = cache
    
    _, slicep, _ = slice_orbitals(Nf, Np, Nv)

    # CONVERT X INTO (A,M)
    X_to_AM!(A, M, X, Nf, Np, Nv)

    # COMPUTE Q(M)
    remove_trace!(QM, M)

    # COMPUTE d = Ω × ([A,ΓM] + Q(M)) × Ω'
    _commutator!(tmp2, A, ΓM, tmp1)
    @views tmp2pp = tmp2[slicep,slicep]
    @. tmp2pp += tmpQM
    _mul!(tmpd, Ω, tmp2, Ω', tmp1)

    # COMPUTE GM = Ω' × (J(d) + Q_xc(DM)) × Ω


    # COMPUTE Z = [FM,A] + GM
    _commutator!(tmpZ, FM, A, tmp1)
    mul!(tmp1, tmpZ, ΓM)
    @. tmpZ += GM

    # COMPUTE   X = 0.5 × ([FM,A]ΓM + FM[A,ΓM]) + GM×ΓM + FM×Q(M) 
    _commutator!(tmp2, A, ΓM, tmp1)
    mul!(tmpX, FM, tmp2)
    @. tmpX += tmp1
    @. tmpX *= 1/2
    mul!(tmp1, GM, ΓM)
    @. tmpX += tmp1
    @views tmp1_p = tmp1[:,slicep]
    @views FM_p   = FM[:,slicep]
    @views tmpX_p = tmpX[:,slicep]
    mul!(tmp1_v,FM_v,tmpQM)
    @. tmpX_v += tmp1_v

    # OUTPUT
    tmp1 .= tmpX .- tmpX'
    @views Zpp = Z[slicep, slicep]
    remove_trace!(Zpp)
    AM_to_X!(vecX, tmp1, Zpp, Nf, Np, Nv)
    vecX
end

function solve_linear_system!(cache::CacheQuadratic, method::Quadratic)

    @unpack vecX, vecB, A, M, Nf, Np, Nv = cache

    L(X::AbstractVector) = linearsystem(cache, X; method.opts...)

    sol = linsolve(L, vecB)

    vecX .= sol[1]

    X_to_AM!(A, M, vecX, Nf, Np, Nv)
end


#####################################################################
#                     UPDATE DENSITY MATRIX
#####################################################################

function update_density!(cache::CacheQuadratic)
    @unpack Ω, Γdiag, M, A, tmp1 = cache
    # UPDATE Ω
    mul!(tmp1, Ω, exp(A))
    Ω .= tmp1
    # UPDATE Λ 
    block(Γdiag)[2] .+= M
end

#####################################################################
#          CONVERSION BETWEEN MATRIX AND VECTOR REPRESENTATIONS
#####################################################################


function X_to_AM!(A::AbstractMatrix, M::AbstractMatrix, X::AbstractVector, Nf::Int, Np::Int, Nv::Int) 
    slicef, slicep, slicev = slice_orbitals(Nf, Np, Nv)
    slicevf, slicevp, slicepf, slicepp = vec_slice_orbitals(Nf, Np, Nv)

    _copy_vec_to_mat!(A, slicev, slicef, X, slicevf)
    _copy_vec_to_mat!(A, slicev, slicep, X, slicevp)
    _copy_vec_to_mat!(A, slicep, slicef, X, slicepf)
    _copy_vec_to_mat!(M, axes(M,1), axes(M,2), X, slicepp)

    LinearAlgebra._copy_adjtrans!(A, slicef, slicev, A, slicev, slicef, antiadjoint)
    LinearAlgebra._copy_adjtrans!(A, slicep, slicev, A, slicev, slicep, antiadjoint)
    LinearAlgebra._copy_adjtrans!(A, slicef, slicep, A, slicep, slicef, antiadjoint)
    nothing
end


function AM_to_X!(X::AbstractVector, A::AbstractMatrix, M::AbstractMatrix, Nf::Int, Np::Int, Nv::Int)
    slicef, slicep, slicev = slice_orbitals(Nf, Np, Nv)
    slicevf, slicevp, slicepf, slicepp = vec_slice_orbitals(Nf, Np, Nv)

    _copy_mat_to_vec!(X, slicevf, A, slicev, slicef)
    _copy_mat_to_vec!(X, slicevp, A, slicev, slicep,)
    _copy_mat_to_vec!(X, slicepf, A, slicep, slicef)
    _copy_mat_to_vec!(X, slicepp, M, axes(M,1), axes(M,2),)
    nothing
end












#=
function solve_system!(solver::KohnShamSolver, method::Quadratic)
    
    @unpack Nf, Np, Nv,
            ΓM, Ω, A, M, B, FA, FM  = cache

    # COMPUTE THE RIGHT-HAND SIDE OF THE SYSTEM
    
    # ???
    setup_one_body_hamiltonian!(FAO, discretization)
    # ???
    
    _mul!(FM, Ω', FA, Ω, tmp1)
    _commutator!(B, FM, ΓM, tmp1)
    _copy_mat_to_vec!(vecB, axes(vecB,1), B, axes(B,1), axes(B,2))

    L = let Nf = Nf, Np = Np, Nv = Nv,
            ΓM = ΓM, Ω = Ω, A = A, M = M, 
            tmp1 = tmp1, tmp2 = tmp2, tmpQM = tmpQM, tmpX = tmpX

            _, slicep, _ = slice_orbitals(Nf, Np, Nv)

            # SETUP THE LINEAR MAP
            function linear_map!(L::AbstractVector, X::AbstractVector)

                # CONVERT X INTO (A,M)
                X_to_AM!(A, M, X, Nf, Np, Nv)

                # COMPUTE Q(M)
                remove_trace!(tmpQM, M)

                # COMPUTE d = Ω × ([A,ΓM] + Q(M)) × Ω'
                _commutator!(tmp2, A, ΓM, tmp1)
                mul!(tmpX, FM, tmp2)
                @views tmp2pp = tmp2[slicep,slicep]
                @. tmp2pp += tmpQM
                _mul!(tmpd, Ω, tmp2, Ω', tmp1)

                # COMPUTE GM = Ω' × (J(d) + Q_xc(DM)) × Ω

                # COMPUTE Z = [FM,A] + GM
                _commutator!(tmpZ, FM, A, tmp1)
                mul!(tmp1, tmpZ, ΓM)
                @. tmpZ += GM

                # COMPUTE   X = 0.5 × ([FM,A]ΓM + FM[A,ΓM]) + GM×ΓM + FM×Q(M)          
                @. tmpX += tmp1
                @. tmpX *= 1/2
                mul!(tmp1, GM, ΓM)
                @. tmpX += tmp1
                @views tmp1_p = tmp1[:,slicep]
                @views FM_p   = FM[:,slicep]
                @views tmpX_p = tmpX[:,slicep]
                mul!(tmp1_v,FM_v,tmpQM)
                @. tmpX_v += tmp1_v

                # OUTPUT
                tmp1 .= tmpX .- tmpX'

                @views Zpp = Z[slicep, slicep]
                AM_to_X!(L, tmp1, Zpp, Nf, Np, Nv)
                nothing
            end

            dim_system = Nv * Nf + Nv * Np + Np * Nf + Np * Np
            LinearOperator( Float64, dim_system, dim_system, false, false,
                                linear_map!, nothing, nothing)
        end

    # SOLVE THE LINEAR SYSTEM
end
=#