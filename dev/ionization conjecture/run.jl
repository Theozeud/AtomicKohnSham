include("../../benchmarktools/atoms/setup.jl")

using AtomicKohnSham
import AtomicKohnSham: solve!
using GLMakie
# PARAMETERS
Rmax = 1000
Nmesh = 50
ordermax =10
lₕ = 1
Z = 3
N = 3
Q = N-Z

X = LinRange(0.0001,1000,10000)
Δx = X[2]-X[1]

# DISCRETIZATION
mesh = expmesh(0,Rmax,Nmesh;s=1.2)
basis = P1IntLegendreBasis(mesh; ordermax = ordermax)
dis = KSEDiscretization(lₕ, basis, mesh, 1, 10)

# MODEL
model = KSEModel(; z = Z, N = N)

# SOLVER
solver = KSESolver(model, dis, ODA(0.4);
            scftol = 1e-10,
            maxiter = 300,
            degen_tol = 1e-10,
            verbose = 0)
# RESOLUTION
solve!(solver)
sol = KSESolution(solver, "")

# POTENTIAL
coeffhart = solver.discretization.cache.tmp_C
Vnuc = @. - Z/X
Vhart0 = evaluate(basis, coeffhart, X)./ X .+ N/Rmax
V = zeros(length(X),lₕ+1)
for l ∈ 0:lₕ
    @views Vl = V[:,l+1]
    @. Vl = Vhart0 + l*(l+1)/(2*X^2) + Vnuc
end


# ORBITAL


# DENSITY
ρX = eval_density(sol, X)

# TOTAL DENSITY FROM R
function cumint_fd2(u::Vector{<:Real}, Δx::Real)
    Nx = length(u)
    I = similar(u)
    @inbounds begin
        I[1] = (Δx/2) * u[1]
        for i in 2:Nx
            I[i] = I[i-1] + (Δx/2) * (u[i-1] + u[i])
        end
    end
    return I
end

g = @. 4π * X^2 * ρX
QX = cumint_fd2(g, Δx)


#=
# PLOT STOPPING CRITERIA
p = plot_stopping_criteria([sol])

# PLOT DENSITY
ρ = plot_density([sol], X)

# PLOT EIGENVECTOR
E = plot_eigenvector(sol,X)
=#
