#####################################################################
#                          LSDA Cache
#####################################################################

mutable struct LSDACache
    A               # Matrix of Qᵢ'Qⱼ'
    M₀              # Matrix of QᵢQⱼ
    M₋₁             # Matrix of 1/x QᵢQⱼ
    M₋₂             # Matrix of 1/x^2 QᵢQⱼ
    F               # Tensor of 1/x QᵢQⱼQₖ
    B               # Matrix of 4πρQᵢ
    C               # Matrix for V(ρ) - solution of the Gauss Electrostatic law
    Cᵨ              # Number equal to total charge
    Kin             # Matrix VF of Kinetic 
    Coulomb         # Matrix VF of Coulomb
    Hfix            # Part of Hamilotnian not needing to be recomputed (Kinetic + Colombial)         
    Hartree         # Matrix VF of hartree 
    Vxc             # Matrix VF of UP   Echange-correlation potential
    Energy          # Total Energy
    Energy_kin      # Kinetic Energy
end

mutable struct LSDA_tmp_Cache
    tmp_H           # Store Hₗ
    tmp_D           # Store ρ
    tmp_Dstar       # Store ρ
    tmp_U           # Store the eigenvector
    tmp_ϵ           # Store the eigenvalue
    tmp_n           # Store the occupation number
    tmp_MV          # Store the contraction F:C  
    tmp_index_sort
end

function create_cache_lsda(lₕ, Nₕ, T)
    A           = zeros(T, Nₕ, Nₕ) 
    M₀          = zeros(T, Nₕ, Nₕ)
    M₋₁         = zeros(T, Nₕ, Nₕ)
    M₋₂         = zeros(T, Nₕ, Nₕ)
    F           = zeros(T, Nₕ, Nₕ, Nₕ)
    B           = zeros(T, Nₕ)
    C           = zeros(T, Nₕ)
    Cᵨ          = zero(T)
    Kin         = zeros(T, lₕ+1, Nₕ, Nₕ)
    Coulomb     = zeros(T, lₕ+1, Nₕ, Nₕ)
    Hfix        = zeros(T, lₕ+1, Nₕ, Nₕ)
    Hartree     = zeros(T, Nₕ, Nₕ)
    Vxc         = zeros(T, Nₕ, Nₕ, 2)
    Energy      = zero(T)
    Energy_kin  = zero(T)

    # Initialization of array for temporary stockage of computations
    tmp_H           = zeros(T, lₕ+1, Nₕ, Nₕ, 2)
    tmp_D           = zeros(T, Nₕ, Nₕ, 2)
    tmp_Dstar       = zeros(T, Nₕ, Nₕ, 2)
    tmp_U           = zeros(T, lₕ+1, Nₕ, Nₕ, 2)
    tmp_MV          = zeros(T, Nₕ, Nₕ)
    tmp_ϵ           = zeros(T, lₕ+1, Nₕ, 2)
    tmp_n           = zeros(T, lₕ+1, Nₕ, 2)   
    tmp_index_sort  = zeros(Int, 2*Nₕ*(lₕ+1))
    LSDACache(A, M₀, M₋₁, M₋₂, F, B, C, Cᵨ, Kin, Coulomb, Hfix, Hartree, Vxc, Energy, Energy_kin),  
    LSDA_tmp_Cache(tmp_H, tmp_D, tmp_Dstar, tmp_U,  tmp_ϵ, tmp_n, tmp_MV, tmp_index_sort)
end


#####################################################################
#                          LSDA Discretization
#####################################################################

struct LSDADiscretization{T} <: KohnShamDiscretization
    lₕ::Int
    Nₕ::Int
    basis::Basis
    mesh::Mesh{T}
    Rmin::T
    Rmax::T
    elT::Type
    cache::LSDACache
    tmp_cache::LSDA_tmp_Cache
    function LSDADiscretization(lₕ::Int, basis::Basis, mesh::Mesh)
        elT = try
             bottom_type(basis)
        catch
            eltype(basis)
        end
        Nₕ = length(basis)
        new{eltype(mesh)}(lₕ, Nₕ, basis, mesh, first(mesh), last(mesh), elT, create_cache_lsda(lₕ, Nₕ, elT)...)
    end
end

#####################################################################
#                          Init Cache
#####################################################################

function init_cache!(discretization::LSDADiscretization, model::AbstractDFTModel, hartree)

    @unpack lₕ, basis  = discretization
    @unpack A, M₀, M₋₁, M₋₂, F, Kin, Coulomb, Hfix = discretization.cache

    # Creation of the base matrices
    A   .= stiffness_matrix(basis)
    M₀  .= mass_matrix(basis)
    M₋₁ .= weight_mass_matrix(basis, -1)

    # Creation of the fix part of the hamiltonian   
    if lₕ ≠ 0
        M₋₂ .= weight_mass_matrix(basis, -2)
    end
    kinetic_matrix!(discretization)
    coulomb_matrix!(discretization, model)
    @. Hfix = Kin + Coulomb

    # Creation of the 3-index tensor F if there is the hartree term
    if !iszero(hartree)
        F .= weight_mass_3tensor(basis, Monomial(-1))
    end

    nothing
end


#####################################################################
#                          Init for Solver
#####################################################################

function init_density_matrix(kd::LSDADiscretization)
    D = zeros(kd.elT,  kd.Nₕ, kd.Nₕ, 2)
    D[:,:,2] .= Matrix(I,kd.Nₕ, kd.Nₕ)
    D[:,:,1] .= Matrix(I,kd.Nₕ, kd.Nₕ)
    return D
end
init_coeffs_discretization(kd::LSDADiscretization) = zeros(kd.elT, kd.lₕ+1, kd.Nₕ, kd.Nₕ, 2)
init_energy(kd::LSDADiscretization)                = zeros(kd.elT, kd.lₕ+1, kd.Nₕ, 2)
init_occupation(kd::LSDADiscretization)            = zeros(kd.elT, kd.lₕ+1, kd.Nₕ, 2)

#####################################################################
#               Find Orbital : Solve the eigen problems
#####################################################################

function find_orbital!(discretization::LSDADiscretization, solver::KohnShamSolver)

    @unpack lₕ = discretization
    @unpack M₀, Hfix, Hartree, Vxc = discretization.cache
    @unpack tmp_H, tmp_U, tmp_ϵ = discretization.tmp_cache
    @unpack Dprev = solver
    @unpack hartree = solver.opts

    # STEP 1 : Compute Hartree term 
    if !iszero(hartree)
        hartree_matrix!(discretization, Dprev)
        @. Hartree = hartree * Hartree
    end

    # STEP 2 : Compute Exchange Correlation term
    if isthereExchangeCorrelation(solver.model)
        exchange_corr_matrix!(discretization, solver.model, Dprev)
    end

    # STEP 3 : Solve the generalized eigenvalue problem for each section l and spin up and down
    for σ ∈ 1:2
        for l ∈ 0:lₕ
            # building the hamiltonian of the lᵗʰ section
            @. tmp_H[l+1,:,:,σ] = Hfix[l+1,:,:] + Vxc[:,:,σ] + Hartree
            # solving
            tmp_ϵ[l+1,:,σ], tmp_U[l+1,:,:,σ] = solve_generalized_eigenvalue_problem(tmp_H[l+1,:,:,σ], M₀)
        end
    end
end

#####################################################################
#               Normamisation of eigenvector
#####################################################################

function normalization!(discretization::LSDADiscretization)
    @unpack M₀ = discretization.cache
    @unpack tmp_U, tmp_n = discretization.tmp_cache
    @unpack lₕ, Nₕ  = discretization
    @inbounds for σ ∈ 1:2
        @inbounds for l ∈ 1:lₕ+1
            @inbounds for k ∈ 1:Nₕ
                if !iszero(tmp_n[l,k,σ])
                    normalization = sqrt(tmp_U[l,:,k,σ]' * M₀ * tmp_U[l,:,k,σ])
                    tmp_U[l,:,k,σ] .= tmp_U[l,:,k,σ] .* 1/normalization
                end
            end
        end
    end
    nothing
end

#####################################################################
#                          Kinetic Matrix
#####################################################################

function kinetic_matrix!(discretization::LSDADiscretization)
    @unpack A, M₋₂, Kin = discretization.cache
    for l ∈ 0:discretization.lₕ
        @. Kin[l+1,:,:] =  1/2 * (A + l*(l+1)*M₋₂)
    end 
    nothing
end

#####################################################################
#                          Coulomb Matrix
#####################################################################

function coulomb_matrix!(discretization::LSDADiscretization, model)
    @unpack M₋₁, Coulomb = discretization.cache
    for l ∈ 0:discretization.lₕ
        Coulomb[l+1,:,:] .= - model.z .* M₋₁
    end 
    nothing
end

#####################################################################
#                          Hartree Matrix
#####################################################################

function hartree_matrix!(discretization::LSDADiscretization, D)
    @unpack A, M₀, F, B, C, Hartree = discretization.cache
    @unpack tmp_MV = discretization.tmp_cache
    @unpack basis, Rmax = discretization
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    DD = DUP .+ DDOWN 
    @tensor B[m] = DD[i,j] * F[i,j,m]
    C .= A\B
    @tensor newCᵨ = DD[i,j] * M₀[i,j]
    @tensor tmp_MV[i,j] = C[k] * F[i,j,k]
    @. Hartree = tmp_MV + newCᵨ/Rmax * M₀
    discretization.cache.Cᵨ = newCᵨ
    nothing
end

#####################################################################
#                   Exchange Correlation Matrix
#####################################################################

function exchange_corr_matrix!(discretization::LSDADiscretization, model, D)
    exchange_corr_matrixUP!(discretization, model, D)
    exchange_corr_matrixDOWN!(discretization, model, D)
    nothing
end

function exchange_corr_matrixUP!(discretization::LSDADiscretization, model, D)
    @unpack Vxc = discretization.cache
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
    weightUP(x) = vxcUP(model.exc, ρUP(x), ρDOWN(x))
    @views VxcUP = Vxc[:,:,1]
    fill_weight_mass_matrix!(discretization.basis, weightUP, VxcUP)
    nothing
end

function exchange_corr_matrixDOWN!(discretization::LSDADiscretization, model, D)
    @unpack Vxc = discretization.cache
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
    weightDOWN(x) = vxcDOWN(model.exc, ρUP(x), ρDOWN(x))
    @views VxcDOWN = Vxc[:,:,2]
    fill_weight_mass_matrix!(discretization.basis, weightDOWN, VxcDOWN)
    nothing
end

#####################################################################
#                             Energy
#####################################################################

function compute_energy!(discretization::LSDADiscretization, solver::KohnShamSolver)
    compute_kinetic_energy!(discretization,solver)
    compute_coulomb_energy!(discretization,solver)
    compute_hartree_energy!(discretization,solver)
    if isthereExchangeCorrelation(solver.model)
        compute_exchangecorrelation_energy!(discretization,solver)
    end
    if isLSDA(solver.model)
        compute_kinetic_correlation_energy!(discretization,solver)
    end
    compute_total_energy!(discretization,solver)
end

function compute_total_energy!(discretization::LSDADiscretization, solver::KohnShamSolver)
    @unpack Rmax = discretization
    @unpack B, C, Cᵨ, Vxc = discretization.cache
    @unpack n, ϵ, D = solver
    @tensor energy = n[l,n,σ] * ϵ[l,n,σ] 
    if isthereExchangeCorrelation(solver.model)
        @tensor energy_correction = Vxc[i,j,σ] * D[i,j,σ]
        solver.energy = energy - discretization.elT(0.5) * (dot(B,C) + Cᵨ^2/Rmax) + solver.energy_exc - energy_correction
    else
        solver.energy = energy - discretization.elT(0.5) * (dot(B,C) + Cᵨ^2/Rmax)
    end
    nothing
end

function compute_kinetic_energy!(discretization::LSDADiscretization, solver::KohnShamSolver)
    @unpack A, M₋₂ = discretization.cache
    @unpack U, n = solver
    @unpack lₕ, Nₕ, elT  = discretization
    solver.energy_kin = zero(solver.energy_kin)
    @inbounds for σ ∈ 1:2
        @inbounds for l ∈ 1:lₕ+1   
            @inbounds for k ∈ 1:Nₕ
                if !iszero(n[l,k,σ])
                    solver.energy_kin += n[l,k,σ] * elT(0.5) * U[l,:,k,σ]' * (A + l*(l+1)*M₋₂) * U[l,:,k,σ]
                end
            end
        end
    end
    nothing
end

function compute_coulomb_energy!(discretization::LSDADiscretization, solver::KohnShamSolver)
    @unpack M₋₁ = discretization.cache
    @unpack U, n = solver
    @unpack lₕ, Nₕ  = discretization
    solver.energy_cou = zero(solver.energy_cou)
    @inbounds for σ ∈ 1:2
        @inbounds for l ∈ 1:lₕ+1   
            @inbounds for k ∈ 1:Nₕ
                if !iszero(n[l,k,σ])
                    solver.energy_cou -= solver.model.z * n[l,k,σ] * U[l,:,k,σ]' * M₋₁ * U[l,:,k,σ]
                end
            end
        end
    end
    nothing
end

function compute_hartree_energy!(discretization::LSDADiscretization, solver::KohnShamSolver)
    @unpack Rmax, elT = discretization
    @unpack B, C, Cᵨ = discretization.cache
    solver.energy_har = elT(0.5) * (dot(B,C) + Cᵨ^2/Rmax)
    nothing
end

function compute_exchangecorrelation_energy!(discretization::LSDADiscretization, solver::KohnShamSolver)
    @unpack D = solver
    @unpack Rmax = discretization
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
    f(x,p) = exc(solver.model.exc, ρUP(x), ρDOWN(x)) * x^2
    prob = IntegralProblem(f, (zero(Rmax),Rmax))
    solver.energy_exc = 4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
    nothing
end

function compute_kinetic_correlation_energy!(discretization::LSDADiscretization, solver::KohnShamSolver)
    @unpack D = solver
    @unpack Rmax = discretization
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    ρUP(x) = compute_densityUP(discretization, DUP, x)
    ρDOWN(x) = compute_densityDOWN(discretization, DDOWN, x)
    ρ(x) = ρDOWN(x) + ρUP(x)
    ξ(x) = (ρUP(x) - ρDOWN(x))/ρ(x)
    rs(x) = (3/(4π * ρ(x)))^(1/3)
    tc(x,p) =  -4 * ec(solver.model.exc, ρUP(x), ρDOWN(x)) * ρ(x) * x^2 + 3 * x^2 * ( ρUP(x)* vcUP(solver.model.exc, ρUP(x), ρDOWN(x))+ ρDOWN(x) * vcDOWN(solver.model.exc, ρUP(x), ρDOWN(x))) #
    prob = IntegralProblem(tc, (zero(Rmax),Rmax))
    solver.energy_kincor = 4π * solve(prob, QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u
    nothing
end

#####################################################################
#                             Density
#####################################################################

function density_matrix!(discretization::LSDADiscretization)
    @unpack tmp_Dstar, tmp_U, tmp_n = discretization.tmp_cache
    @unpack lₕ, Nₕ  = discretization
    @inbounds for σ ∈ 1:2
        @inbounds for l ∈ 1:lₕ+1 
            @inbounds for k ∈ 1:Nₕ
                if !iszero(tmp_n[l,k,σ])
                    @inbounds for i ∈ 1:Nₕ
                        val = tmp_n[l,k,σ] * tmp_U[l,i,k,σ] 
                        @inbounds @simd for j ∈ 1:i
                            tmp_Dstar[i,j,σ] += val * tmp_U[l,j,k,σ]
                        end
                    end
                end
            end
        end
    end
    @inbounds for i in 1:Nₕ
        @inbounds @simd for j in 1:i-1
            tmp_Dstar[j,i,:] = tmp_Dstar[i,j,:]
        end
    end
    nothing
end

function compute_density(discretization::LSDADiscretization, D, x)
    @unpack basis = discretization
    newT = promote_type(eltype(basis), typeof(x))
    val = zero(newT)
    eval_basis = zeros(newT, basis.size)
    @inbounds for i ∈ eachindex(basis)
        eval_basis[i] = basis(i,x)
    end
    @views DUP = D[:,:,1]   
    @views DDOWN = D[:,:,2]
    val = (eval_basis)' * (DUP + DDOWN) * eval_basis
    return val* 1/4π * 1/(x^2)
end

function compute_densityUP(discretization::LSDADiscretization, DUP, x)
    @unpack basis = discretization
    newT = promote_type(eltype(basis), typeof(x))
    val = zero(newT)
    eval_basis = zeros(newT, basis.size)
    @inbounds for i ∈ eachindex(basis)
        eval_basis[i] = basis(i,x)
    end
    val = (eval_basis)' * DUP * eval_basis
    return val* 1/4π * 1/(x^2)

end

function compute_densityDOWN(discretization::LSDADiscretization, DDOWN, x)
    @unpack basis = discretization
    newT = promote_type(eltype(basis), typeof(x))
    val = zero(newT)
    eval_basis = zeros(newT, basis.size)
    @inbounds for i ∈ eachindex(basis)
        eval_basis[i] = basis(i,x)
    end
    val = (eval_basis)' * DDOWN * eval_basis
    return val* 1/4π * 1/(x^2)
end