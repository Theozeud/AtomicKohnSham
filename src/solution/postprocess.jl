# ===================================================================
#                        EVALUATION OF EIGENVECTORS
# ===================================================================
"""
    eval_orbital(sol, n, l, X; h10=false)
    eval_orbital(sol, n, l, σ, X; h10=false)
    eval_orbital(sol, idx, X; h10=false)

Evaluate a Kohn–Sham orbital on the FEM basis at radial points `X`.

This function evaluates the radial part of a Kohn–Sham orbital associated with
the quantum numbers `(n, l)` (and spin `σ` if applicable) using the FEM basis
stored in the solution context.

By default, the returned orbital corresponds to the physical radial wavefunction
`u(r) / r`. If `h10=true`, the function instead returns the FEM representation
`u(r)` in the space `H¹₀`, without division by `r`.

# Arguments
- `sol::KSESolution`: Converged Kohn–Sham solution.
- `n::Int`: Principal quantum number.
- `l::Int`: Orbital angular momentum quantum number.
- `σ::Int`: Spin index (required for spin-polarized calculations).
- `idx::String`: Shell label (e.g. `"1s"`, `"2p↑"`), parsed internally.
- `X::AbstractVector`: Radial evaluation points.

# Keyword Arguments
- `h10::Bool`: If `true`, return the FEM `H¹₀` representation; otherwise return
  the physical radial orbital.

# Returns
A vector containing the orbital values evaluated at points `X`.
"""
function eval_orbital(sol::KSESolution, n::Int, l::Int, X::AbstractVector{<:Real};
                     h10::Bool = false)
    @unpack model, basis = context
    @assert 0 ≤ l ≤ n-1 "Wrong number quantum. You should have 0 ≤ l ≤ n-1."
    @assert model.nspin == 1 "The discretization is spin-polarized. Please give a spin σ."
    φh10X = evaluate(basis, sol.orbitals[:, n - l, l + 1], X)
    if h10
        return φh10X
    else
        return φh10X ./ X
    end
end

function eval_orbital(sol::KSESolution, n::Int, l::Int, σ::Int, X::AbstractVector{<:Real};
                     h10::Bool = false)
    @assert 0 ≤ l ≤ n-1 "Wrong number quantum. You should have 0 ≤ l ≤ n-1."
    φh10X = evaluate(sol.context.basis, sol.orbitals[:, n - l, l + 1, σ], X)
    if h10
        return φh10X
    else
        return φh10X ./ X
    end
end

function eval_orbital(sol::KSESolution, idx::String, X::AbstractVector{<:Real};
                     h10::Bool = false)
    qn = parse_shell(idx)
    eval_orbital(sol, qn..., X; h10=h10)
end

# ===================================================================
#                       EVALUATION OF THE DENSITY
# ===================================================================
"""
    eval_density(sol, X)

Evaluate the (spin-summed) electron density on radial points `X`.

For spin-unpolarized calculations (`nspin == 1`), this returns `ρ(X)` computed
from the density matrix stored in `sol`. For spin-polarized calculations
(`nspin == 2`), this returns the total density `ρ↑(X) + ρ↓(X)`.

# Arguments
- `sol::KSESolution`: Converged Kohn–Sham solution.
- `X::AbstractVector`: Radial evaluation points.

# Returns
A vector containing the density values evaluated at `X`.
"""
function eval_density(sol::KSESolution{T},
                      X::AbstractVector{TX}) where {T<:Real, TX <: Real}
    @unpack model, basis = sol.context
    newT = promote_type(T,TX)
    if model.nspin == 1
        ρ = zeros(newT, length(X))
        eval_density!(ρ, basis, sol.D, X)
        return ρ
    else
        ρup = zeros(newT, length(X))
        @views DUP = sol.D[:, :, 1]
        eval_density!(ρup, basis, DUP, X)
        ρdown = zeros(newT, length(X))
        @views DDOWN = sol.D[:, :, 2]
        eval_density!(ρdown, basis, DDOWN, X)
        @. ρup += ρdown
        return ρup
    end
end

"""
    eval_density(sol, X, σ)

Evaluate the spin density `ρσ` on radial points `X`.

This method is available for spin-polarized calculations and returns the density
associated with spin channel `σ`.

# Arguments
- `sol::KSESolution`: Converged Kohn–Sham solution.
- `X::AbstractVector`: Radial evaluation points.
- `σ::Int`: Spin index (1-based).

# Returns
A vector containing `ρσ(X)` evaluated at `X`.
"""
function eval_density(sol::KSESolution{T}, X::AbstractVector{TX},
                      σ::Int) where {T<:Real, TX <: Real}
    @unpack model, basis = sol.context
    @assert 1 ≤ σ ≤ model.nspin
    @views Dσ = sol.D[:, :, σ]
    newT = promote_type(T,TX)
    ρσ = zeros(newT, length(X))
    eval_density!(ρσ, basis, Dσ, X)
    ρσ
end


"""
    eval_density!(ρ, basis, D, X)

Evaluate the radial electron density at points `X` from a density matrix `D`,
storing the result in `ρ`.

This in-place routine avoids allocations and is intended for internal use
during quadrature and energy evaluations.
"""
function eval_density!(ρ::AbstractVector{T},
                      basis::FEMBasis,
                      D::AbstractMatrix{<:Real},
                      X::AbstractVector{<:Real}) where T
    buf1 = zeros(T, basis.max_nb_poly_cells)
    buf2 = zeros(T, basis.max_nb_poly_cells)
    @inbounds for k in eachindex(X)
        xk = X[k]
        localisation_xk = findindex(basis.mesh, xk)
        Ik = basis.cells_to_indices[localisation_xk]
        @views eval_basis = buf2[Ik]
        evaluate!(eval_basis, basis, Ik, xk)
        @views tv = buf1[Ik]
        @views Dk = D[Ik, Ik]
        mul!(tv, Dk, eval_basis)
        ρ[k] = 1/(4π*xk^2) * dot(eval_basis, tv)
    end
    ρ
end


# ===================================================================
#                    EVALUATION OF THE POTENTIALS
# ===================================================================
function eval_hartree_pot(sol::KSESolution{T},
                         X::AbstractVector{TX}) where {T<:Real, TX <: Real}

    @unpack model, basis = sol.context
    Vhart = evaluate(basis, W, X)./ X .+ model.N/last(basis.mesh)
    return Vhart
end
