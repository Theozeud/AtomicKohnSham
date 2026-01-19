using AtomicKohnSham
using LinearAlgebra
using GenericLinearAlgebra
using UnPack
import AtomicKohnSham: FunWeight, fill_mass_vector!

# ==========================================================================================
#                               Hydrogenoid FEM utilities
# ==========================================================================================
"""
Defines a hydrogenoid radial FEM problem.
"""
struct HydrogenoidProblem{T<:Real,Fmesh,Fbasis}
    Z::T
    l::Int
    eigs_idx::Vector{Int}
    Rmax::T
    Nmesh::Int
    mesh_type::Fmesh
    mesh_opts::NamedTuple
    basis_type::Fbasis
    basis_opts::NamedTuple
    name::String
end

Base.eltype(::HydrogenoidProblem{T}) where T = T

function HydrogenoidProblem(; T::Type{<:Real},
                            Z::Real,
                            l::Int,
                            Rmax::Real,
                            Nmesh::Int,
                            mesh_type,
                            mesh_opts::NamedTuple = NamedTuple(),
                            basis_type,
                            basis_opts::NamedTuple = NamedTuple(),
                            eigs_idx = collect(l+1:l+3),
                            name::AbstractString = "")
    return HydrogenoidProblem{T,typeof(mesh_type),typeof(basis_type)}(
        T(Z), l, eigs_idx, T(Rmax), Nmesh,
        mesh_type, mesh_opts,
        basis_type, basis_opts,
        String(name),
    )
end

"""
    HydrogenoidProblem(prob; kwargs...)

Copy-constructor: returns a new HydrogenoidProblem by copying `prob`
and replacing only the fields specified in the keyword arguments.
"""
function HydrogenoidProblem(prob::HydrogenoidProblem;
                            T        = eltype(prob),
                            Z        = prob.Z,
                            l        = prob.l,
                            Rmax     = prob.Rmax,
                            Nmesh    = prob.Nmesh,
                            mesh_type  = prob.mesh_type,
                            mesh_opts  = prob.mesh_opts,
                            basis_type = prob.basis_type,
                            basis_opts = prob.basis_opts,
                            eigs_idx = prob.eigs_idx,
                            name     = prob.name)
    return HydrogenoidProblem(; T=T, Z=Z, l=l,
                                Rmax=Rmax, Nmesh=Nmesh,
                                mesh_type=mesh_type, mesh_opts=mesh_opts,
                                basis_type=basis_type, basis_opts=basis_opts,
                                eigs_idx=eigs_idx, name=name)
end


"""
    HydrogenoidSolution

Container for:
- selected numerical eigenvalues/eigenvectors
- theoretical eigenvalues/eigenvectors
- error estimates Δλ and ΔU
"""
struct HydrogenoidSolution{T<:Real}
    problem::HydrogenoidProblem
    eigs_idx::Vector{Int}
    λ::Vector{T}
    U::Matrix{T}
    λtheo::Vector{T}
    Utheo
    Δλ::Vector{T}
    ΔU::Vector{T}
end

function HydrogenoidSolution(problem::HydrogenoidProblem,
                             λ::AbstractVector{<:Real},
                             U::AbstractMatrix{<:Real},
                             A::AbstractMatrix{<:Real},
                             M0::AbstractMatrix{<:Real},
                             basis)
    @unpack eigs_idx,l = problem
    λ_sel = λ[eigs_idx .- l]
    λtheo = theoretical_eigenvalue(problem, eigs_idx)
    Δλ = abs.(λ_sel .- λtheo)

    U_sel = U[:, eigs_idx.- l]
    Utheo = theoretical_eigenvector(problem, eigs_idx)
    V = zeros(length(basis))
    C = zero(V)
    invM0 = inv(M0)
    ΔU = zeros(length(eigs_idx))
    for i ∈ eachindex(eigs_idx)
        fun = Utheo[i]
        weight = FunWeight(fun;is_vectorized=false,is_inplace=false)
        fill_mass_vector!(basis,V;weight=weight)
        mul!(C, invM0,V)
        h1norm_pos = (U_sel[:,i]-C)' * (A+M0) * (U_sel[:,i]-C)
        h1norm_neg = (-U_sel[:,i]-C)' * (A+M0) * (-U_sel[:,i]-C)
        if h1norm_neg < h1norm_pos
            @. U_sel[:,i] *= -1
        end
        ΔU[i] = min(h1norm_pos,h1norm_neg)
    end

    return HydrogenoidSolution{eltype(problem)}(problem, eigs_idx,
                                  λ_sel, U_sel,
                                  λtheo, Utheo,
                                  Δλ, ΔU)
end



# ==========================================================================================
#                                   Assembly & solver
# ==========================================================================================
"""
    L2normalization(U, M0)

Normalizes each column of matrix `U` in the L² radial norm induced by
the mass matrix `M0`.
"""
function L2normalization(U::AbstractMatrix, M0::AbstractMatrix)
    @assert size(U,1) == size(M0,1) == size(M0,2)
    U_norm = similar(U)
    M0uk = zero(U[:,1])
    for k in axes(U, 2)
        uk = view(U, :, k)
        mul!(M0uk,M0,uk)
        norm2 = dot(uk, M0uk)
        normalization = sqrt(norm2)
        U_norm[:, k] .= uk ./ normalization
    end
    return U_norm
end


"""
    assemble_hydro_operators(problem)

Builds the FEM Hamiltonian H and mass matrix M0 for the radial
hydrogenoid operator:

    H u = -1/2 u'' + l(l+1)/(2 r^2) u - Z/r * u

Returns (H, M0, basis).
"""
function assemble_hydro_operators(problem::HydrogenoidProblem)
    @unpack Z, l, Rmax, Nmesh, mesh_type, mesh_opts, basis_type, basis_opts = problem
    T = eltype(problem)
    mesh  = mesh_type(zero(T), Rmax, Nmesh; T=T, mesh_opts...)
    basis = basis_type(mesh, T; basis_opts...)
    A   = Symmetric(stiffness_matrix(basis))
    M0  = Symmetric(mass_matrix(basis))
    M_1 = Symmetric(mass_matrix(basis, -1))
    if l == 0
        H = T(0.5) * A - Z * M_1
        return H, A, M0, basis
    else
        M_2 = Symmetric(mass_matrix(basis, -2))
        H = T(0.5) * (A + l*(l+1) * M_2) - Z * M_1
        return H, A, M0, basis
    end
end


"""
    eigen_hydro(problem)

Solves the generalized eigenvalue problem:

    H u = λ M0 u

Performs L² radial normalization and returns a HydrogenoidSolution.
"""
function eigen_hydro(problem::HydrogenoidProblem)
    H, A, M0, basis = assemble_hydro_operators(problem)
    Λ, U = eigen(H, M0)
    U_norm = L2normalization(U, M0)
    #return Λ, U_norm
    return HydrogenoidSolution(problem, Λ, U_norm, A, M0, basis)
end


# ==========================================================================================
#                           Theoretical eigenvectors and eigenvalue
# ==========================================================================================
"""
    theoretical_eigenvector(n, Z, l, T = Float64)

Return a function `u(r)` giving the theoretical hydrogenoid radial
solution u_{n,l}(r) for the equation

    -1/2 u''(r) + l(l+1)/(2r^2) u(r) - Z/r * u(r) = λ u(r),

with eigenvalue λ_n = -Z^2 / (2 n^2).

Here u(r) is the "folded" radial function: u(r) = r * R_{n,l}(r),
normalized in the radial L² sense (∫ |u(r)|² dr = 1).
"""
function theoretical_eigenvector(n::Int, Z::Real, l::Int, T = Float64)
    @assert n > l ≥ 0 "Quantum numbers must satisfy n > l ≥ 0."

    # Degree of associated Laguerre polynomial: k = n - l - 1
    k = n - l - 1
    α = 2l + 1

    # Build coefficients of L_k^{(α)}(x) = ∑_{m=0}^k c_m x^m
    # c_m = (-1)^m * binomial(k + α, k - m) / m!
    c = zeros(T, k + 1)
    for m in 0:k
        c[m+1] = (-1)^m * binomial(k + α, k - m) / factorial(big(m))
    end

    # Normalization constant for R_{nl}(r):
    # R_{nl}(r) = N_{nl} e^{-ρ/2} ρ^l L_k^{(2l+1)}(ρ),  ρ = 2Zr/n
    # N_{nl} = 2 Z^{3/2} / n^2 * sqrt( (n-l-1)! / (n+l)! )
    N_radial = 2 * Z^(3/2) / (n^2) * sqrt(factorial(big(n - l - 1)) / factorial(big(n + l)))

    return function unl(r::Real)
        ρ = 2 * Z * r / n
        # Evaluate L_k^{(α)}(ρ) by Horner's method
        L = zero(T)
        for m in k:-1:0
            L = L * ρ + c[m+1]
        end
        # R_{nl}(r)
        R = N_radial * exp(-ρ / 2) * ρ^l * L
        # u(r) = r * R_{nl}(r)
        u = r * R
        return u
    end
end

"""
    theoretical_eigenvector(problem, idx)

Return a vector of functions [u_i(r)] corresponding to the theoretical
eigenvectors for the indices in `idx` (assuming hydrogenic model).
"""
function theoretical_eigenvector(problem::HydrogenoidProblem{T},
                                 idx::AbstractVector{Int}) where T
    @unpack Z, l = problem
    return [theoretical_eigenvector(n, Z, l, T) for n in idx]
end


function theoretical_eigenvalue(problem::HydrogenoidProblem{T},
                                idx::AbstractVector{Int}) where T
    @unpack Z, l = problem
    return [-Z^2/(2*(n+l)^2) for n in idx]
end
