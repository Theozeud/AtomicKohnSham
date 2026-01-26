"""
    KSESolver(model::KSEModel,
              discretization::KSEDiscretization,
              alg::SCFAlgorithm;
              maxiter::Int = 100,
              logconfig::LogConfig = LogConfig(),
              verbose::Int = 0) -> KSESolver

Create a solver for the radial Kohn–Sham extended (KSE) equations, using a given
`model`, `discretization`, and SCF algorithm `alg`.

A `KSESolver` holds all state required to run an iterative SCF procedure:
the discretization and model, the chosen SCF algorithm and its workspace (cache),
the current energy components, and logging/monitoring data.

Run the SCF loop with `solve!(solver)` and build a result object with
`KSESolution(solver)`.

# Arguments
- `model::KSEModel`: Physical model.
- `discretization::KSEDiscretization`: Numerical discretization.
- `alg::SCFAlgorithm`: SCF strategy.

# Keyword arguments
- `maxiter::Int = 100`: Maximum number of SCF iterations.
- `logconfig::LogConfig = LogConfig()`: Logging configuration.
- `verbose::Int = 0`: Verbosity level
  (`0`: silent, `1`: iterations, `2`: + basic metrics, `3`: + method-specific details).

# Fields
- `niter::Int`: Number of completed SCF iterations.
- `stopping_criteria::T`: Current value of the convergence metric.
- `discretization`: Discretization object.
- `model`: Model object (physics / functionals / parameters).
- `alg`: SCF algorithm instance (options/parameters).
- `algcache`: Algorithm workspace (mutable cache used across iterations).
- `energies::Energies{T}`: Current energy decomposition.
- `logbook`: Logging state and recorded history.
- `maxiter::Int`: Maximum number of iterations.
- `verbose::UInt8`: Stored verbosity level.
"""
mutable struct KSESolver{
    D<:KSEDiscretization,
    M<:KSEModel,
    A<:SCFAlgorithm,
    C<:SCFCache,
    L<:LogBook,
    T<:Real,
    B
}
    niter::Int
    stopping_criteria::T
    discretization::D
    model::M
    alg::A
    algcache::C
    energies::Energies{T}
    logbook::L
    callback::B
    maxiter::Int
    verbose::UInt8

    function KSESolver(model::M, discretization::D, alg::A; maxiter::Int = 100,
                      verbose::Int = 0, callback::B = nothing) where {D<:KSEDiscretization,
                        M<:KSEModel, A<:SCFAlgorithm, B}
        # Data type of numbers
        T = eltype(discretization)
        # Initial numbers
        niter = 0
        stopping_criteria = zero(T)
        energies = Energies(T)
        # Init cache of the discretization
        init_cache!(discretization, model)
        # Create the cache of the algorithm
        algcache = create_cache_alg(alg, discretization, model)
        # Create the logbook
        logbook = LogBook(T)
        return new{D, M, A, typeof(algcache), typeof(logbook), T, B}(
            niter,
            stopping_criteria,
            discretization,
            model,
            alg,
            algcache,
            energies,
            logbook,
            callback,
            Int(maxiter),
            UInt8(verbose))
    end
end
