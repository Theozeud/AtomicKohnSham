"""
    KSESolver(model::KSEModel,
              discretization::KSEDiscretization,
              alg::SCFAlgorithm;
              maxiter::Int = 100,
              callback = nothing) -> KSESolver

Create a solver for the radial Kohn–Sham extended (KSE) equations, using a given
`model`, `discretization`, and SCF algorithm `alg`.

A `KSESolver` holds all state required to run an iterative SCF procedure:
the discretization and model, the chosen SCF algorithm and its workspace (cache),
the current energy components, and logging data.

Run the SCF loop with `solve!(solver)` and build a result object with
`KSESolution(solver)`. For per-iteration diagnostics, pass a callback (e.g.
[`LogFileCallback`](@ref), wrapped in a [`CallbackSet`](@ref)) rather than a
verbosity level — it can write to `stdout` just as well as to a file.

# Arguments
- `model::KSEModel`: Physical model.
- `discretization::KSEDiscretization`: Numerical discretization.
- `alg::SCFAlgorithm`: SCF strategy.

# Keyword arguments
- `maxiter::Int = 100`: Maximum number of SCF iterations.
- `callback = nothing`: Callback (or [`CallbackSet`](@ref)) invoked once per
  SCF iteration.

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

    function KSESolver(model::M, discretization::D, alg::A; maxiter::Int = 100,
                      callback::B = nothing) where {D<:KSEDiscretization,
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
            Int(maxiter))
    end
end
