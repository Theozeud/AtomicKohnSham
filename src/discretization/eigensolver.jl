abstract type EigenSolver end

"""
    FullEigenSolver()

Diagonalize the whole generalized eigenproblem and keep only the lowest `nₕ`
eigenpairs. Always correct, but wastes O(Nₕ³) work when `nₕ ≪ Nₕ`.
"""
struct FullEigenSolver <: EigenSolver end

"""
    PartialEigenSolver()

Ask LAPACK (`syevr` via `eigen(A, 1:nₕ)`) directly for the lowest `nₕ`
eigenpairs, without diagonalizing the rest of the spectrum. Default choice:
cheaper than [`FullEigenSolver`](@ref) whenever `nₕ < Nₕ`.
"""
struct PartialEigenSolver <: EigenSolver end

solve_eigenproblem(::FullEigenSolver, A::Symmetric, nₕ::Int) = eigen(A)
solve_eigenproblem(::PartialEigenSolver, A::Symmetric, nₕ::Int) = eigen(A, 1:nₕ)
