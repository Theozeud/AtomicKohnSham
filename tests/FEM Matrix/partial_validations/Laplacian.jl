using KohnShamResolution
using LinearAlgebra
using Plots
using UnPack
using GenericLinearAlgebra
using DoubleFloats

using IterativeSolvers

# Structure Problem
struct LaplacianProblem
    name
    T
    typebasis
    optsbasis
    typemesh
    optsmesh
    Rmax
    Nmesh

    LaplacianProblem(;name, T, typebasis, optsbasis, typemesh, optsmesh, Rmax, Nmesh) =
        new(name, T, typebasis, optsbasis, typemesh, optsmesh, Rmax, Nmesh)

    LaplacianProblem(prob;  name = prob.name, T = prob.T, typebasis = prob.typebasis, 
                            optsbasis = prob.optsbasis, typemesh = prob.typemesh, 
                            optsmesh = prob.optsmesh, Rmax = prob.Rmax, Nmesh = prob.Nmesh) =
        new(name, T, typebasis, optsbasis, typemesh, optsmesh, Rmax, Nmesh)
end

# Structure Solution
struct LaplacianSolution
    prob
    λ
    Δλ
    u
    function LaplacianSolution(prob, λ, u)
        λ_theo = theoretical_eigenvalue.(eachindex(λ), prob.Rmax, prob.T) 
        Δλ = abs.(λ .- λ_theo)
        new(prob, λ, Δλ, u)
    end
end

# Structure to compute convergence curve
struct LaplacianConvergenceNmesh
    probs           # Set of problems
    vecNmesh        # Set of Nmesh used
    ΔΛ              # Dict of errors on eigenvalues : for each problem,
                    # there is a vector of errors depending on Nmesh
    nums            # Number of eigenvalue and eigenvector used 
end

# Theoretical Eigenvalue and Eigenvectors
theoretical_eigenvector(n, Rmax, T) = x-> sqrt(T(2)/T(Rmax)) * sin(n*T(π)*x/T(Rmax)) 
theoretical_eigenvalue(n, Rmax, T)  = (T(π)/T(Rmax))^2 * n^2

# Compute eigenvectors and eigenvalues
function eigen_(problems)
    sols = LaplacianSolution[]
    for problem ∈ problems
        @unpack T, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis = problem
        mesh = typemesh(zero(T), Rmax, Nmesh; T = T, optsmesh...)
        basis = typebasis(mesh, T; optsbasis...)
        M₀  = Symmetric(mass_matrix(basis))
        A   = Symmetric(stiffness_matrix(basis))
        u, λ = eigen(A, M₀)
        push!(sols, LaplacianSolution(problem, λ, u))
    end
    sols
end

# Compute eigen values (usefull for Double64 for which eigen is not defined)
function eigvals_(problems)
    sols = LaplacianSolution[]
    for problem ∈ problems
        @unpack T, Rmax, Nmesh, typemesh, optsmesh, typebasis, optsbasis = problem
        mesh = typemesh(zero(T), Rmax, Nmesh; T = T, optsmesh...)
        basis = typebasis(mesh, T; optsbasis...)
        M₀  = mass_matrix(basis)
        A   = Symmetric(stiffness_matrix(basis))
        λ = eigvals(inv(M₀) * A)
        #λ  = lobpcg(A, M₀, false, 10; tol = 1e-9).λ
        push!(sols, LaplacianSolution(problem, λ, nothing))
    end
    sols
end



# Plot Eigenvalue
function plot_eigenvalue(sols)
    plt = plot( size= (1500, 700), margin = 0.5Plots.cm,                
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)
    mindex = 0
    RMAX   = 0
    for sol ∈ sols
        plot!(plt, sol.λ, label = sol.prob.name, lw = 5)
        if length(sol.λ) > mindex
            mindex = length(sol.λ)
        end
        if sol.prob.Rmax > RMAX
            RMAX = sol.prob.Rmax
        end
    end
    plot!(plt, theoretical_eigenvalue.(1:mindex, RMAX, Float64), label = "Theoretical", lw = 4, ls = :dash, lc = :black)
    return plt
end

##############################################################################
# Plot Eigenvector

function plot_eigenvector(sols, n; resol = 5000)
    plt = plot( size = (1500, 700), margin = 1Plots.cm,                
                legendfontsize  = 14,  
                titlefontsize   = 14,
                guidefontsize   = 14,
                tickfontsize    = 14)
    xlabel!(plt, "r")
    title!(plt, string(n)*"-th eigenvector")
    RMAX = 0
    for sol ∈ sols
        eign = build_on_basis(basis, sol.u[:,n])
        eign_n = eign / sqrt(integrate(eign*eign, zero(T), T(Rmax)))
        X = LinRange(O, sol.prob.Rmax, resol)
        plot!(plt, X, (sign(theoretical_eigenvector(n, Rmax, sol.prob.T)(1) * eign_n(1)) * eign_n).(X), label = sol.prob.name, lw = 4)
        if sol.prob.Rmax > RMAX
            RMAX = sol.prob.Rmax
        end
    end
    X = LinRange(O, RMAX, resol)
    plot!(plt, X, theoretical_eigenvector(n, Rmax, sol.prob.T).(X), label = "Theoretical n = $n", lw = 4)
    plt
end

##################################################################################
# Compute Error

function convergenceNmesh(vecNmesh::AbstractVector, problems; nums = [1])
    ΔΛ = Dict()
    for problem ∈ problems
        @unpack T, name = problem
        println(name)
        Δλ = zeros(T, length(vecNmesh), length(nums))
        @inbounds for i ∈ eachindex(vecNmesh)
            newprob = LaplacianProblem(problem; Nmesh = vecNmesh[i])
            @time "Nmesh = $(vecNmesh[i])" sol = eigvals_([newprob])
            Δλ[i,:] .= sol[1].Δλ[nums]
        end
        ΔΛ[name] = Δλ
    end
    LaplacianConvergenceNmesh(problems, vecNmesh, ΔΛ, nums)
end

function convergence_plot_Nmesh(sols::LaplacianConvergenceNmesh)
    plt = plot( size = (1000,800), margin = 0.5Plots.cm, legend = :bottomleft , xaxis=:log, yaxis=:log,
                legendfontsize  = 14,  
                titlefontsize   = 18,
                guidefontsize   = 14,
                tickfontsize    = 14)
    xlabel!(plt, "Nmesh")
    ylabel!(plt, "Error")
    title!(plt, "Convergence Plot")
    for prob ∈ sols.probs
        for num ∈ sols.nums
            plot!(plt, sols.vecNmesh, sols.ΔΛ[prob.name][:,num], lw = 4, label = prob.name*"-$num", markershape = :x, markersize = 10)
        end
    end
    plt
end

##############################################
prob = LaplacianProblem(name = "problem", 
                        T = Float64,
                        typebasis = P1IntLegendreGenerator,
                        optsbasis = (ordermax = 2,),
                        typemesh = linmesh,
                        optsmesh = (;),
                        Rmax = 100,
                        Nmesh = 100)

sol = eigvals_([prob])
plot_eigenvalue(sol)

sol = convergenceNmesh(2 .^(4:8), [prob])

convergence_plot_Nmesh(sol)