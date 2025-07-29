#--------------------------------------------------------------------
#                  STRUCTURE OF THE SOLUTION
#--------------------------------------------------------------------
"""
    KSESolution(solver::KSESolver, name::String) -> KSESolution

Constructs a solution object from a fully converged (or stopped) Kohn–Sham solver.
This structure stores the result of a self-consistent field (SCF) computation,
including metadata such as convergence status, final energies, iteration logs,
and algorithm-specific data.

# Arguments
- `solver::KSESolver`: The SCF solver after a call to `solve!`.
- `name::String`: Name or label for the solution (can be empty if unused).

# Fields
- `success::String`: Convergence status — either `"SUCCESS"` if convergence occurred before reaching `maxiter`,
                                            or   `"MAXITERS"` if the iteration limit was reached.
- `solveropts::SolverOptions`: Options that were used in the solver (e.g., tolerance, integration method).
- `niter::Int`: Number of SCF iterations performed.
- `stopping_criteria::T`: Final convergence metric value (e.g., norm of residual density).
- `energies::Dict{Symbol, T}`: Energies computed (total, kinetic, Hartree, etc.).
- `datas::SCFSolution`: Struct with detailed solution data, whose fields depend on the SCF algorithm used.
- `log::LogBook`: Object that records iteration history and diagnostics.
- `name::String`: Name of the solution (can be used for identification, labeling results, etc.).
"""
struct KSESolution{ optsType <: SolverOptions,
                    T <: Real,
                    solutionType <: SCFSolution,
                    logbookType <: LogBook,
                    D <: KSEDiscretization}

    success::String                     # Print the final state of the solver
                                        #   Can be : SUCCESS or MAXITERS

    solveropts::optsType                # Option of the solver used

    niter::Int                          # Number of iterations
    stopping_criteria::T                # Final stopping criteria

    energies::Dict{Symbol, T}           # Energies

    datas::solutionType                 # All datas depending on the scf algorithm used

    log::logbookType                    # LogBook of all recorded quantities

    name::String                        # Name of the solution

    discretization::D                   # Discretization used to compute the solution

    function KSESolution(solver::KSESolver, name::String)

        # FLAG ON THE SUCCESS (OR NOT) OF THE MINIMIZATION
        success = solver.niter == solver.opts.maxiter ? "MAXITERS" : "SUCCESS"

        # DATAS
        datas = makesolution(solver.cache, solver.alg, solver)
        discretization = solver.discretization
        new{typeof(solver.opts),
            typeof(solver.stopping_criteria),
            typeof(datas),
            typeof(solver.logbook),
            typeof(discretization)}(success,
                                    solver.opts,
                                    solver.niter,
                                    solver.stopping_criteria,
                                    solver.energies,
                                    datas,
                                    solver.logbook,
                                    name,
                                    discretization)
    end
end


function Base.getproperty(sol::KSESolution, s::Symbol)
    if s ∈ fieldnames(KSESolution)
        getfield(sol, s)
    elseif s ∈ propertynames(getfield(sol, :datas))
        getfield(getfield(sol, :datas), s)
    else
        throw(ErrorException("type KSESolution has no field $(s)"))
    end
end


#--------------------------------------------------------------------
#                  DISPLAY A SUMMARY OF THE SOLUTION
#--------------------------------------------------------------------
const L_QUANTUM_LABELS = ("s", "p", "d", "f", "g", "h", "i")
const SPIN_LABELS = ("↑", "↓")

function Base.show(io::IO, sol::KSESolution)
    printstyled(io, "Name : " * (sol.name) * "\n"; bold = true)
    printstyled(io, "Sucess = "; bold = true)
    if sol.success == "SUCCESS"
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :green)
    else
        printstyled(io, string(sol.success)*"\n"; bold = true, color = :red)
    end
    printstyled(io, "niter = "; bold = true)
    println(io, string(sol.niter))
    printstyled(io, "Stopping criteria = "; bold = true)
    println(io, string(sol.stopping_criteria ))
    printstyled(io, "All Energies :\n"; bold = true, color = :green)
    for s ∈ keys(sol.energies)
        printstyled(io, "            $(s) = $(sol.energies[s]) \n"; bold = true, color = :green)
    end
    printstyled(io, "Occupation number = \n"; bold = true, color = :blue)
    for i ∈ eachindex(sol.occupation_number)
        occupation_number = sol.occupation_number[i]
        printstyled(io, "            $(occupation_number[1]) : ($(occupation_number[2]),$(occupation_number[3])) \n"; bold = true, color = :blue)
    end
end


#--------------------------------------------------------------------
#                  POST-PROCESSING COMPUTATIONS
#--------------------------------------------------------------------

# COMPUTE EIGENVECTORS
function eigenvector(
    sol::KSESolution,
    n::Int,
    l::Int,
    X::AbstractVector{<:Real})
    @assert 0 ≤ l ≤ n-1 "Wrong number quantum. You should have 0 ≤ l ≤ n-1."
    @assert sol.discretization.n_spin == 1 "The discretization is spin-polarized. Please give a spin σ."
    evaluate(sol.discretization.basis, sol.orbitals[:, n-l,l+1], X)
end

function eigenvector(sol::KSESolution, n::Int, l::Int, σ::Int, X::AbstractVector{<:Real})
    @assert 0 ≤ l ≤ n-1 "Wrong number quantum. You should have 0 ≤ l ≤ n-1."
    evaluate(sol.discretization.basis, sol.orbitals[:, n-l,l+1,σ], X)
end



# COMPUTE DENSITY
function eval_density(
    sol::KSESolution,
    X::AbstractVector{<:Real})
    if sol.discretization.n_spin == 1
        eval_density(sol.discretization, sol.density_coeffs, X)
    else
        @views DUP = sol.density_coeffs[:,:,1]
        ρUP = eval_density(sol.discretization, DUP, X)
        @views DDOWN = sol.density_coeffs[:,:,1]
        ρDOWN = eval_density(sol.discretization, DDOWN, X)
        ρUP .+ ρDOWN
    end
end


function eval_density(
    sol::KSESolution,
    X::AbstractVector{<:Real},
    σ::Int)
    @assert 1≤ σ ≤ sol.discretization.n_spin
    @views Dσ = sol.density_coeffs[:,:,σ]
    eval_density(sol.discretization, Dσ, X)
end


# TOTAL CHARGE OF THE SYSTEM
function total_charge(sol::KSESolution)
    @unpack Rmax = sol.problem.discretization
    f(x,p) = density(sol, x) * x^2
    prob = IntegralProblem(f, (zero(Rmax),Rmax))
    return 4π * solve(prob, QuadGKJL(); reltol = 1e-13, abstol = 1e-13).u
end
