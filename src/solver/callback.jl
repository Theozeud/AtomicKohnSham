"""
    CallbackSet

Composite callback container.

Groups multiple callbacks into a single object. At each SCF iteration, all
registered callbacks are executed sequentially.
"""
struct CallbackSet{C}
    callbacks::C
end

CallbackSet() = CallbackSet(())

"""
Apply all callbacks to the current solver state.
"""
callback!(cb::CallbackSet, solver::KSESolver) =
    foreach(c -> callback!(c, solver), cb.callbacks)

callback!(::Nothing, solver::KSESolver) = nothing
