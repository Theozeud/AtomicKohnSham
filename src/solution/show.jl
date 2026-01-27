# Nombre de chiffres significatifs "raisonnable" dérivé de eps(T)
# (Float64 -> ~16, Float32 -> ~7, etc.)
_sig_digits(::Type{T}) where {T} = T <: AbstractFloat ? max(1, floor(Int, -log10(eps(one(T)))) + 1) : 6

# Largeur "safe" de l'exposant (pour couvrir e-..., e+...)
_exp_digits(::Type{T}) where {T} = T <: AbstractFloat ? length(string(abs(exponent(floatmax(T))))) : 3

# Format scientifique sans Printf :
# ex: "-1.2345e-123"
function _sci_str(x::T; sig::Int=_sig_digits(T)) where {T}
    # Gestion des cas simples
    if x == 0
        return "0.0e+0"
    end
    if !isfinite(x)
        return string(x)  # "Inf", "-Inf", "NaN"
    end

    sgn = x < 0 ? "-" : ""
    ax  = abs(x)

    # Exposant décimal
    e = floor(Int, log10(ax))
    # Mantisse dans [1,10)
    m = ax / (T(10)^e)

    # Arrondi à sig chiffres significatifs
    # (on arrondit sur la mantisse)
    if sig > 1
        scale = T(10)^(sig - 1)
        m = round(m * scale) / scale
        # Si l'arrondi fait passer à 10.0, on renormalise
        if m ≥ 10
            m /= 10
            e += 1
        end
    else
        m = round(m)
        if m ≥ 10
            m /= 10
            e += 1
        end
    end

    # On force un point décimal pour la stabilité d'affichage
    ms = string(m)
    if !occursin(".", ms)
        ms *= ".0"
    end

    esign = e ≥ 0 ? "+" : "-"
    return sgn * ms * "e" * esign * string(abs(e))
end

# Largeur de colonne "ε = ..." déterminée automatiquement (marge incluse)
function _eps_col_width(::Type{T}; extra=4) where {T}
    sig = _sig_digits(T)
    ed  = _exp_digits(T)

    # worst-case approx: "-d." + (sig-1 digits) + "e" + sign + expdigits
    # => 1(sign) + 1(digit) + 1(dot) + (sig-1) + 1(e) + 1(exp sign) + ed
    base = 1 + 1 + 1 + max(sig-1, 0) + 1 + 1 + ed
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

    #=
    printstyled(io, "All Energies :\n"; bold = true, color = :green)
    for s in keys(sol.energies)
        printstyled(io, "            $(s) = $(sol.energies[s]) \n"; bold = true, color = :green)
    end
    =#
    printstyled(io, "Occupation number = \n"; bold = true, color = :blue)

    sig = _sig_digits(T)
    eps_col_width = _eps_col_width(T; extra=4)

    for occ in sol.occupation_number
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
