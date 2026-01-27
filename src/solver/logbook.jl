"""
    LogBook{T}

Lightweight container for monitoring SCF convergence.

Stores the history of stopping criteria and total energies during the SCF
iterations. This structure is intended for diagnostics and convergence analysis.
"""
struct LogBook{T <: Real}
    stopping_criteria::Vector{T}
    Etot::Vector{T}
    function LogBook(::Type{T}) where T
        stopping_criteria   = T[]
        Etot                = T[]
        new{T}(stopping_criteria, Etot)
    end
end
