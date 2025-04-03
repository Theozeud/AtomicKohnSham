#####################################################################
#                  STRUCTURE OF THE SOLUTION
#####################################################################

const L_QUANTUM_LABELS = ("s", "p", "d", "f", "g", "h", "i")
const SPIN_LABELS = ("↑", "↓")

struct KohnShamSolution{problemType <: DFTProblem, 
                        optsType <: SolverOptions, 
                        T <: Real, 
                        solutionType <: SCFSolution,
                        logbookType <: LogBook}

    problem::problemType                # Problem solved

    success::String                     # Print the final state of the solver
                                        #        Can be : SUCCESS or MAXITERS

    solveropts::optsType                # Option of the solver used
    
    niter::Int                          # Number of iterations
    stopping_criteria::T                # Final stopping criteria
              
    energies::Dict{Symbol, T}           # Energies

    datas::solutionType                 # All datas depending on the method used

    log::logbookType                    # LogBook of all recorded quantities

    name::String                        # Name of the solution

    function KohnShamSolution(solver::KohnShamSolver, name::String, datas::SCFSolution)

        # CREATION OF A PROBLEM STRUCTURE TO STORE THE PROBLEM SOLVED
        problem = DFTProblem(solver.model, solver.discretization, solver.method)

        # FLAG ON THE SUCCESS (OR NOT) OF THE MINIMIZATION
        success = solver.niter == solver.opts.maxiter ? "MAXITERS" : "SUCCESS"

        new{typeof(problem),
            typeof(solver.opts),
            typeof(solver.stopping_criteria),
            typeof(datas),
            typeof(solver.logbook)}(problem, 
                                    success, 
                                    solver.opts, 
                                    solver.niter, 
                                    solver.stopping_criteria, 
                                    solver.energies,
                                    datas,
                                    solver.logbook, 
                                    name)
    end
end


function Base.getproperty(sol::KohnShamSolution, s::Symbol)
    if s ∈ fieldnames(KohnShamSolution)
        getfield(sol, s)
    elseif s ∈ propertynames(getfield(sol, :datas))
        getfield(getfield(sol, :datas), s)
    else
        throw(ErrorException("type KohnShamSolution has no field $(s)"))
    end
end


#####################################################################
#                  DISPLAY A SUMMARY OF THE SOLUTION
#####################################################################

function Base.show(io::IO, sol::KohnShamSolution)
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
        display_occupation_number(io, sol.problem.discretization, sol.occupation_number[i])
    end
end

function display_occupation_number(io::IO, ::LDADiscretization, occupation_number)
    printstyled(io, "            $(occupation_number[1]) : ($(occupation_number[2]),$(occupation_number[3])) \n"; bold = true, color = :blue)
end

function display_occupation_number(io::IO, ::LSDADiscretization, occupation_number)
    printstyled(io, "            $(occupation_number[1]) : ($(occupation_number[2]),$(occupation_number[3])) \n"; bold = true, color = :blue)
end

#####################################################################
#                  POST-PROCESSING COMPUTATIONS
#####################################################################

# COMPUTE EIGENVECTORS

function eigenvector(sol::KohnShamSolution, n::Int, l::Int, σ::Int, x)
    eigenvector(sol.problem.discretization, sol, n, l, σ, x)
end

function eigenvector(discretization::LDADiscretization, sol::KohnShamSolution, n::Int, l::Int, σ::Int, x)
    @assert 0 ≤ l ≤ n-1
    tmp = discretization.basis(sol.orbitals[l+1,:, n-l], x)
    if iszero(x) && tmp ≈ zero(tmp)
        return zero(tmp)
    else
        return  1/sqrt(4π * x) * tmp 
    end
end

function eigenvector(discretization::LSDADiscretization, sol::KohnShamSolution, n::Int, l::Int, σ::Int, x)
    @assert 0 ≤ l ≤ n-1
    tmp = discretization.basis(sol.orbitals[l+1,:, n-l,σ], x)
    if iszero(x) && tmp ≈ zero(tmp)
        return zero(tmp)
    else
        return  1/sqrt(4π * x) * tmp 
    end
end


# COMPUTE DENSITY

function density(sol::KohnShamSolution, x::Real)
    compute_density(sol.problem.discretization, sol.density_coeffs, x)
end

function density(sol::KohnShamSolution, X::AbstractVector)
    [compute_density(sol.problem.discretization, sol.density_coeffs, x) for x ∈ X]
end

# TOTAL CHARGE OF THE SYSTEM
function total_charge(sol::KohnShamSolution)
    @unpack Rmax = sol.problem.discretization
    f(x,p) = density(sol, x) * x^2
    prob = IntegralProblem(f, (zero(Rmax),Rmax))
    return 4π * solve(prob, QuadGKJL(); reltol = 1e-13, abstol = 1e-13).u
end