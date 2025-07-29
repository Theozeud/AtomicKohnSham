using JuMP
using OSQP
using LinearAlgebra

"""
    project_trace_blocks(M_blocks, alphas)

Projette une matrice diagonale par blocs `M_blocks` sur l'ensemble des matrices par blocs
dont la trace totale est nulle et chaque trace partielle est >= alpha.

# Entrées
- `M_blocks` : vecteur de matrices carrées `M_ell` (par exemple, Vector{Matrix{Float64}})
- `alphas` : vecteur des α_ell (même longueur que `M_blocks`)

# Sortie
- `M_projected` : vecteur de matrices projetées (même forme que M_blocks)
"""
function project_trace_blocks(M_blocks::Vector{Matrix{Float64}}, alphas::Vector{Float64})
    L = length(M_blocks)
    d = [size(M, 1) for M in M_blocks]
    tr_M = [tr(M) for M in M_blocks]

    # Setup the QP
    model = Model(OSQP.Optimizer)
    set_silent(model)

    @variable(model, t[1:L]) # projected traces
    @objective(model, Min, sum((t[ℓ] - tr_M[ℓ])^2 / d[ℓ] for ℓ in 1:L))

    @constraint(model, sum(t) == 0)
    @constraint(model, [ℓ=1:L], t[ℓ] >= alphas[ℓ])

    x = zero(alphas[1])
    for i in 1:(L - 1)
        set_start_value(t[i], alphas[i])
        x += alphas[i]
    end
    set_start_value(t[L], -x)

    optimize!(model)

    # Extract optimal traces
    t_proj = value.(t)

    # Project each block to have trace t_proj[ℓ]
    M_proj = [M - (tr(M) - t_proj[ℓ]) / d[ℓ] * I for (ℓ, M) in enumerate(M_blocks)]
    return M_proj
end

function init_t(α::AbstractVector{<:Real})
    t0 = copy(α)
    t0[end] = - sum(α) + α[end]
    t0
end

function check_constrains(t::AbstractVector{<:Real}, α::AbstractVector{<:Real})
    all(t .≥ α)
end

# Exemple : 3 blocs
M1 = randn(2, 2);
M1 = (M1 + M1')/2
M2 = randn(3, 3);
M2 = (M2 + M2')/2
M3 = randn(1, 1)

M_blocks = [M1, M2, M3]
alphas = [-1.0, 0.0, 0.6]

M_proj = project_trace_blocks(M_blocks, alphas)

# Vérification
println("Traces projetées : ", [tr(M) for M in M_proj])
println("Somme des traces : ", sum(tr.(M_proj)))  # ≈ 0
println("Contraintes alpha respectées ? ", all(tr.(M_proj) .>= alphas))
