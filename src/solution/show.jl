using Printf

"Number of significant digits worth printing for a float type `T`."
_sig_digits(::Type{T}) where {T} = max(1, floor(Int, -log10(eps(one(T)))) + 1)

"Format `x` in scientific notation with `sig` significant digits."
function _sci_str(x::T; sig::Int = _sig_digits(T)) where {T}
    if x == 0
        return "0.0e+0"
    end
    if !isfinite(x)
        return string(x)
    end
    prec = max(sig - 1, 0)
    s = @sprintf("%.*e", prec, x)
    return s
end

"Column width needed to print any value of type `T` in scientific notation without wrapping."
function _eps_col_width(::Type{T}; extra=4) where {T}
    sig = _sig_digits(T)
    smax = _sci_str(floatmax(T); sig=sig)
    base = max(length(smax), length("-" * smax))
    return length("ε = ") + base + extra
end

"""
Print the "Occupation number" block (shell label, orbital energy, occupation)
shared by `Base.show(io, sol)` and `write_report`.
"""
function _print_occupied(io::IO, sol::KSESolution{T}; styled::Bool = true) where {T}
    if styled
        printstyled(io, "Occupation number = \n"; bold = true, color = :blue)
    else
        println(io, "Occupation number = ")
    end

    sig = _sig_digits(T)
    eps_col_width = _eps_col_width(T; extra=2)

    for (orb, ϵ, n) in sol.occupied
        eps_part = "ε = " * (T <: AbstractFloat ? _sci_str(T(ϵ); sig=sig) : string(ϵ))
        eps_part = rpad(eps_part, eps_col_width)

        line = "            $(orb) : $(eps_part)     n = $(n)\n"
        if styled
            printstyled(io, line; bold = true, color = :blue)
        else
            print(io, line)
        end
    end
end

"""
Print the "Sanity checks" block (energy balance, electron count, density
normalization, orbital orthonormality, virial ratio) shared by
`Base.show(io, sol)` and `write_report`. See [`SanityChecks`](@ref) for what
each value means and why a nonzero error indicates a real problem.
"""
function _print_sanity(io::IO, sol::KSESolution{T}; styled::Bool = true) where {T}
    if styled
        printstyled(io, "Sanity checks = \n"; bold = true, color = :blue)
    else
        println(io, "Sanity checks = ")
    end
    sc = sol.sanity
    _write_kv_block(io, 12, [
        ("energy balance (Etot - ΣE)", _sci_str(sc.energy_balance)),
        ("electron count error (Σn - N)", _sci_str(sc.electron_count_error)),
        ("density norm error (∫4πr²ρdr - N)", _sci_str(sc.density_norm_error)),
        ("orthonormality error (‖UᵀM₀U - I‖∞)", _sci_str(sc.orthonormality_error)),
        ("virial ratio (Ekin / -Etot)", _sci_str(sc.virial_ratio)),
    ])
end

function Base.show(io::IO, sol::KSESolution{T}) where {T}
    printstyled(io, "Name : $(sol.name)\n"; bold = true)
    printstyled(io, "Success = "; bold = true)
    printstyled(io, "$(sol.success)\n"; bold = true,
        color = sol.success == "SUCCESS" ? :green : :red)
    printstyled(io, "niter = "; bold = true)
    println(io, sol.niter)
    printstyled(io, "Stopping criteria = "; bold = true)
    println(io, sol.stopping_criteria)

    _print_occupied(io, sol)
    println(io)
    _print_sanity(io, sol)
end
