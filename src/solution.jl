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
                    logbookType <: LogBook}

    success::String                     # Print the final state of the solver
                                        #        Can be : SUCCESS or MAXITERS

    solveropts::optsType                # Option of the solver used

    exc::Symbol                         # Type of exchange correlation functional
    
    niter::Int                          # Number of iterations
    stopping_criteria::T                # Final stopping criteria
              
    energies::Dict{Symbol, T}           # Energies

    datas::solutionType                 # All datas depending on the scf algorithm used

    log::logbookType                    # LogBook of all recorded quantities

    name::String                        # Name of the solution

    function KSESolution(solver::KSESolver, name::String)

        # FLAG ON THE SUCCESS (OR NOT) OF THE MINIMIZATION
        success = solver.niter == solver.opts.maxiter ? "MAXITERS" : "SUCCESS"

        # DATAS
        datas = makesolution(solver.cache, solver.alg, solver)

        new{typeof(solver.opts),
            typeof(solver.stopping_criteria),
            typeof(datas),
            typeof(solver.logbook)}(success, 
                                    solver.opts, 
                                    typeexc(solver.model),
                                    solver.niter, 
                                    solver.stopping_criteria, 
                                    solver.energies,
                                    datas,
                                    solver.logbook, 
                                    name)
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
    #printstyled(io, "Energy = $(sol.energies[:Etot]) \n"; bold = true, color = :green)
    printstyled(io, "All Energies :\n"; bold = true, color = :green)
    for s ∈ keys(sol.energies)
        printstyled(io, "            $(s) = $(sol.energies[s]) \n"; bold = true, color = :green)
    end
    #printstyled(io, "            $(s) = $(sol.energies[s]) \n"; bold = true, color = :green)
    printstyled(io, "Occupation number = \n"; bold = true, color = :blue)
    for i ∈ eachindex(sol.occupation_number)
        occupation_number = sol.occupation_number[i]
        printstyled(io, "            $(occupation_number[1]) : ($(occupation_number[2]),$(occupation_number[3])) \n"; bold = true, color = :blue)
    end
end


function display_occupation_number(io::IO, ::KSEDiscretization, occupation_number)
    printstyled(io, "            $(occupation_number[1]) : ($(occupation_number[2]),$(occupation_number[3])) \n"; bold = true, color = :blue)
end

#=
function display_occupation_number(io::IO, ::LSDADiscretization, occupation_number)
    printstyled(io, "            $(occupation_number[1]) : ($(occupation_number[2]),$(occupation_number[3])) \n"; bold = true, color = :blue)
end
=#

#--------------------------------------------------------------------
#                  POST-PROCESSING COMPUTATIONS
#--------------------------------------------------------------------

# COMPUTE EIGENVECTORS

function eigenvector(sol::KSESolution, n::Int, l::Int, σ::Int, x)
    eigenvector(sol.problem.discretization, sol, n, l, σ, x)
end

function eigenvector(discretization::KSEDiscretization, sol::KSESolution, n::Int, l::Int, σ::Int, x)
    @assert 0 ≤ l ≤ n-1
    tmp = discretization.basis(sol.orbitals[l+1,:, n-l], x)
    if iszero(x) && tmp ≈ zero(tmp)
        return zero(tmp)
    else
        return  1/sqrt(4π * x) * tmp 
    end
end

#=
function eigenvector(discretization::LSDADiscretization, sol::KSESolution, n::Int, l::Int, σ::Int, x)
    @assert 0 ≤ l ≤ n-1
    tmp = discretization.basis(sol.orbitals[l+1,:, n-l,σ], x)
    if iszero(x) && tmp ≈ zero(tmp)
        return zero(tmp)
    else
        return  1/sqrt(4π * x) * tmp 
    end
end
=#

# COMPUTE DENSITY

function density(sol::KSESolution, x::Real)
    compute_density(sol.problem.discretization, sol.density_coeffs, x)
end

function density(sol::KSESolution, X::AbstractVector)
    [compute_density(sol.problem.discretization, sol.density_coeffs, x) for x ∈ X]
end

# TOTAL CHARGE OF THE SYSTEM
function total_charge(sol::KSESolution)
    @unpack Rmax = sol.problem.discretization
    f(x,p) = density(sol, x) * x^2
    prob = IntegralProblem(f, (zero(Rmax),Rmax))
    return 4π * solve(prob, QuadGKJL(); reltol = 1e-13, abstol = 1e-13).u
end