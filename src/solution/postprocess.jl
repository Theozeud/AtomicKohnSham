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
    @unpack model, basis = sol.context
    @assert 0 ≤ l ≤ n-1 "Wrong number quantum. You should have 0 ≤ l ≤ n-1."
    @assert model.nspin == 1 "The discretization is spin-polarized. Please give a spin σ."
    _eval_orbital_h10(basis, sol.U[:, n - l, l + 1], X, h10)
end

function eval_orbital(sol::KSESolution, n::Int, l::Int, σ::Int, X::AbstractVector{<:Real};
                     h10::Bool = false)
    @unpack basis = sol.context
    @assert 0 ≤ l ≤ n-1 "Wrong number quantum. You should have 0 ≤ l ≤ n-1."
    _eval_orbital_h10(basis, sol.U[:, n - l, l + 1, σ], X, h10)
end

function eval_orbital(sol::KSESolution, idx::String, X::AbstractVector{<:Real};
                     h10::Bool = false)
    qn = parse_shell(idx)
    eval_orbital(sol, qn..., X; h10=h10)
end

"Evaluate the FEM `H¹₀` representation `u(r)`, dividing by `r` unless `h10`."
function _eval_orbital_h10(basis::FEMBasis, coeffs::AbstractVector{<:Real},
                          X::AbstractVector{<:Real}, h10::Bool)
    φh10X = evaluate(basis, coeffs, X)
    return h10 ? φh10X : φh10X ./ X
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
storing the result in `ρ`. Used internally by [`eval_density`](@ref) for
post-hoc evaluation/plotting of a converged solution (the SCF loop itself
evaluates the density through the unrelated
`eval_density!(ρ, discretization, D, X)` method in
`discretization/density.jl`, which is unaffected by this method).

For each point `xₖ ∈ X`, only the FEM basis functions supported on the mesh
cell containing `xₖ` contribute (all others vanish there by compact support),
so the density is evaluated locally:

    ρ(xₖ) = 1/(4π xₖ²) ∑ᵢⱼ∈Iₖ Dᵢⱼ χᵢ(xₖ) χⱼ(xₖ),

where `Iₖ` are the (global) indices of the basis functions supported on that
cell. `buf1`/`buf2` are small scratch buffers reused across points, sized to
the largest number of basis functions supported on any single cell
(`basis.max_nb_poly_cells`) rather than the full basis dimension `Nₕ` — they
must therefore be indexed by *local* position (`1:length(Iₖ)`), not by the
(global, possibly out-of-range) values in `Iₖ` itself.

# Arguments
- `ρ::AbstractVector`: Output vector, overwritten in place, `length(ρ) == length(X)`.
- `basis::FEMBasis`: FEM basis the density matrix `D` is expressed in.
- `D::AbstractMatrix{<:Real}`: One-body reduced density matrix (`Nₕ × Nₕ`).
- `X::AbstractVector{<:Real}`: Radial evaluation points.
"""
function eval_density!(ρ::AbstractVector{T},
                      basis::FEMBasis,
                      D::AbstractMatrix{<:Real},
                      X::AbstractVector{<:Real}) where T
    buf1 = zeros(T, basis.max_nb_poly_cells)
    buf2 = zeros(T, basis.max_nb_poly_cells)
    cache_Pϕx = _cache_Pϕx(basis, first(X))
    @inbounds for k in eachindex(X)
        xk = X[k]
        localisation_xk = findindex(basis.mesh, xk)
        Ik = basis.cells_to_indices[localisation_xk]
        nk = length(Ik)
        @views eval_basis = buf2[1:nk]
        evaluate!(eval_basis, basis, Ik, xk, cache_Pϕx)
        @views tv = buf1[1:nk]
        @views Dk = D[Ik, Ik]
        mul!(tv, Dk, eval_basis)
        ρ[k] = 1/(4π*xk^2) * dot(eval_basis, tv)
    end
    ρ
end


# ===================================================================
#                    EVALUATION OF THE POTENTIALS
# ===================================================================
"""
    eval_hartree(sol::KSESolution, X::AbstractVector{<:Real})

Evaluate the radial Hartree potential at the radial points `X`.

The Hartree potential solves `-1/r (r Vᴴ)'' = 4πρ` with the finite-domain
boundary lifting `θ(r) = r/Rmax`, giving (see the model derivation)

    Vᴴ(r) = W(r) / r + N / Rmax,

where `W` is the FEM representation of the regular part `w` of the Hartree
potential (stored in `sol.W`, solving `𝔸 W = 𝔹[𝔻]`), `N` is the number of
electrons, and `Rmax` is the last point of the radial mesh. As with
[`eval_orbital`](@ref), `X` should not contain `r = 0` (the FEM representation
`W` vanishes there, so `W(r)/r` needs to be evaluated as a limit).

# Arguments
- `sol::KSESolution`: Converged or stopped Kohn–Sham solution.
- `X::AbstractVector{<:Real}`: Radial evaluation points.

# Returns
A vector containing the Hartree potential values evaluated at `X`.
"""
function eval_hartree(sol::KSESolution,
                      X::AbstractVector{TX}) where {TX<:Real}
    @unpack model, basis = sol.context
    WX = evaluate(basis, sol.W, X)
    return WX ./ X .+ model.N/last(basis.mesh)
end


"""
    eval_nuclear(sol::KSESolution, X::AbstractVector{<:Real})

Evaluate the electron-nucleus attraction potential at the radial points `X`.

The nuclear potential is

    Vₙᵤ꜀(r) = -Z / r,

where `Z` is the nuclear charge stored in `sol.context.model`. At `r = 0`,
the value is `-Inf`.

# Arguments
- `sol::KSESolution`: Kohn–Sham solution containing the model parameters.
- `X::AbstractVector{<:Real}`: Radial evaluation points.

# Returns
A vector containing the nuclear potential values evaluated at `X`.
"""
function eval_nuclear(sol::KSESolution{T},
                      X::AbstractVector{TX}) where {T<:Real, TX<:Real}
    Z = sol.context.model.Z
    Tout = float(promote_type(typeof(Z), TX))
    Vnuc = similar(X, Tout)
    @inbounds for i in eachindex(X)
        x = X[i]
        Vnuc[i] = iszero(x) ? -oftype(one(Tout), Inf) : -Z / x
    end
    return Vnuc
end


"""
    eval_kinetic_potential(l::Int, X::AbstractVector{<:Real})
    eval_kinetic_potential(sol::KSESolution, l::Int, X::AbstractVector{<:Real})

Evaluate the centrifugal part of the radial kinetic operator at the radial
points `X`.

For an angular momentum quantum number `l`, the centrifugal potential is

    Vₖᵢₙ,l(r) = l(l + 1) / (2r²).

At `r = 0`, the value is `0` if `l == 0`, and `+Inf` if `l > 0`. The `sol`
argument is accepted (and ignored) only so this function shares its call
signature with the other potential evaluators (see [`eval_effective_potential`](@ref)).

# Arguments
- `l::Int`: Orbital angular momentum quantum number.
- `X::AbstractVector{<:Real}`: Radial evaluation points.
- `sol::KSESolution`: Optional solution argument, included for API consistency.

# Returns
A vector containing the centrifugal kinetic potential values evaluated at `X`.
"""
function eval_kinetic_potential(l::Int,
                                X::AbstractVector{TX}) where {TX<:Real}
    @assert l ≥ 0 "The angular momentum quantum number l must be non-negative."
    c = l * (l + 1) / 2
    Tout = float(promote_type(typeof(c), TX))
    Vkin = similar(X, Tout)
    @inbounds for i in eachindex(X)
        x = X[i]
        if iszero(x)
            Vkin[i] = iszero(c) ? zero(Tout) : oftype(one(Tout), Inf)
        else
            Vkin[i] = c / x^2
        end
    end
    return Vkin
end

function eval_kinetic_potential(::KSESolution,
                                l::Int,
                                X::AbstractVector{<:Real})
    return eval_kinetic_potential(l, X)
end


"""
    eval_vxc(sol::KSESolution, X::AbstractVector{<:Real})
    eval_vxc(sol::KSESolution, X::AbstractVector{<:Real}, σ::Int)

Evaluate the local (LDA/LSDA) exchange–correlation potential `vₓ꜀ = dεₓ꜀/dρ`
at the radial points `X`.

The density is first evaluated at `X` with [`eval_density`](@ref), then the
exchange–correlation functionals stored in `sol.context.model` are evaluated
pointwise on that density, using the same `evaluate_vrho!` routine as during
SCF assembly (`assemble_exc!`). Returns a vector of zeros if the model has no
exchange–correlation (`has_exchcorr(model) == false`).

The first method is for spin-unpolarized calculations (`nspin == 1`); the
second requires a spin index `σ` and is for spin-polarized calculations,
returning `vₓ꜀,σ(ρ↑(r), ρ↓(r))`.

# Arguments
- `sol::KSESolution`: Converged Kohn–Sham solution.
- `X::AbstractVector{<:Real}`: Radial evaluation points.
- `σ::Int`: Spin index (required for spin-polarized calculations).

# Returns
A vector containing the exchange–correlation potential evaluated at `X`.
"""
function eval_vxc(sol::KSESolution{T}, X::AbstractVector{TX}) where {T<:Real, TX<:Real}
    @unpack model = sol.context
    @assert model.nspin == 1 "The discretization is spin-polarized. Please give a spin σ."
    has_exchcorr(model) || return zeros(promote_type(T, TX), length(X))
    ρX = eval_density(sol, X)
    vxc = similar(ρX)
    cache = similar(ρX)
    evaluate_vrho!(model; rho = ρX, vrho = vxc, cache = cache)
    return vxc
end

function eval_vxc(sol::KSESolution{T}, X::AbstractVector{TX},
                  σ::Int) where {T<:Real, TX<:Real}
    @unpack model = sol.context
    @assert 1 ≤ σ ≤ model.nspin
    has_exchcorr(model) || return zeros(promote_type(T, TX), length(X))
    ρup = eval_density(sol, X, 1)
    ρdown = eval_density(sol, X, 2)
    ρ = permutedims(hcat(ρup, ρdown))
    vxc = similar(ρ)
    cache = similar(ρ)
    evaluate_vrho!(model; rho = ρ, vrho = vxc, cache = cache)
    @views return vxc[σ, :]
end


"""
    eval_effective_potential(sol::KSESolution, l::Int, X::AbstractVector{<:Real}; σ::Int=1)

Evaluate the total effective radial potential seen by an orbital of angular
momentum `l` (and spin `σ` if applicable), i.e. the potential part of the
radial Kohn–Sham mean-field Hamiltonian:

    Vₑff,l,σ(r) = Vₙᵤ꜀(r) + Vₖᵢₙ,l(r) + Vᴴ(r) + Vₓ꜀,σ(r),

summing [`eval_nuclear`](@ref), [`eval_kinetic_potential`](@ref),
[`eval_hartree`](@ref) (only if `model.hartree != 0`), and [`eval_vxc`](@ref)
(only if `has_exchcorr(model)`).

# Arguments
- `sol::KSESolution`: Converged Kohn–Sham solution.
- `l::Int`: Orbital angular momentum quantum number.
- `X::AbstractVector{<:Real}`: Radial evaluation points.

# Keyword Arguments
- `σ::Int`: Spin index, used only for spin-polarized calculations.

# Returns
A vector containing the effective potential values evaluated at `X`.
"""
function eval_effective_potential(sol::KSESolution, l::Int, X::AbstractVector{<:Real};
                                  σ::Int = 1)
    @unpack model = sol.context
    Veff = eval_nuclear(sol, X) .+ eval_kinetic_potential(sol, l, X)
    iszero(model.hartree) || (Veff .+= eval_hartree(sol, X))
    if has_exchcorr(model)
        Veff .+= model.nspin == 1 ? eval_vxc(sol, X) : eval_vxc(sol, X, σ)
    end
    return Veff
end
