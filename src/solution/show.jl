using Printf

_sig_digits(::Type{T}) where {T} = max(1, floor(Int, -log10(eps(one(T)))) + 1)

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

function _eps_col_width(::Type{T}; extra=4) where {T}
    sig = _sig_digits(T)
    smax = _sci_str(floatmax(T); sig=sig)
    base = max(length(smax), length("-" * smax))
    return length("ε = ") + base + extra
end

function Base.show(io::IO, sol::KSESolution{T}) where {T}
    printstyled(io, "Name : " * (sol.name) * "\n"; bold = true)
    printstyled(io, "Sucess = "; bold = true)
    if sol.success == "SUCCESS"
        printstyled(io, string(sol.success) * "\n"; bold = true, color = :green)
    else
        printstyled(io, string(sol.success) * "\n"; bold = true, color = :red)
    end
    printstyled(io, "niter = "; bold = true)
    println(io, string(sol.niter))
    printstyled(io, "Stopping criteria = "; bold = true)
    println(io, string(sol.stopping_criteria))

    printstyled(io, "Occupation number = \n"; bold = true, color = :blue)

    sig = _sig_digits(T)
    eps_col_width = _eps_col_width(T; extra=2)

    for occ in sol.occupied
        orb = occ[1]
        ϵ   = occ[2]
        n   = occ[3]

        eps_part = "ε = " * (T <: AbstractFloat ? _sci_str(T(ϵ); sig=sig) : string(ϵ))
        eps_part = rpad(eps_part, eps_col_width)

        printstyled(io,
            "            $(orb) : $(eps_part)     n = $(n)\n";
            bold = true, color = :blue)
    end
end
