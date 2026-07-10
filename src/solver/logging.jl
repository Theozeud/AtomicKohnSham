# ===================================================================
#                       FILE LOGGING CALLBACK
# ===================================================================
"""
    LogFileCallback(io::IO)

Callback that appends one line per SCF iteration to `io`, meant to be
combined into a [`CallbackSet`](@ref) and passed to `groundstate`/`KSESolver`
via the `callback` keyword. This replaces ad-hoc `redirect_stdout` logging:
`io` is owned and closed by the caller, and logging keeps working even if the
SCF loop throws.

# Example
```julia
open("log.txt", "w") do io
    write_log_header(io, model, alg, dis)
    sol = groundstate(model, dis, alg; maxiter = 100,
                       callback = CallbackSet((LogFileCallback(io),)))
end
```
"""
struct LogFileCallback{IOT<:IO}
    io::IOT
end

function callback!(cb::LogFileCallback, solver::KSESolver)
    @unpack io = cb
    @unpack niter, stopping_criteria, energies, algcache = solver
    t = getfield_or_nothing(algcache, :t)
    line = "iter = $(niter)   stop = $(stopping_criteria)   Etot = $(energies.Etot)" *
        "   Ekin = $(energies.Ekin)   Ecou = $(energies.Ecou)   Ehar = $(energies.Ehar)" *
        "   Eexc = $(energies.Eexc)"
    line = t === nothing ? line : line * "   t = $(t)"
    println(io, line)
    flush(io)
    nothing
end

"""
    write_log_header(io::IO, model::KSEModel, alg::SCFAlgorithm, discretization::KSEDiscretization)

Write the parameter block (identical to the one in [`write_report`](@ref))
at the top of a `log.txt` file, before the per-iteration lines appended by
[`LogFileCallback`](@ref).
"""
function write_log_header(io::IO, model::KSEModel, alg::SCFAlgorithm,
                         discretization::KSEDiscretization)
    @unpack lₕ, nₕ, basis, fem_integration_method = discretization
    _write_parameters(io, model, alg, lₕ, nₕ, basis, fem_integration_method)
    println(io)
    println(io, "ITERATIONS")
    println(io, "----------")
    nothing
end
