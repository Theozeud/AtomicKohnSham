const L_QUANTUM_LABELS = ('s','p','d','f','g','h','i')
const SPIN_LABELS = ('↑','↓')   # σ=1 -> ↑, σ=2 -> ↓

# -------------------------
# Helpers
# -------------------------
@inline function l_from_label(c::Char)::Int
    i = findfirst(==(c), L_QUANTUM_LABELS)
    i === nothing && error("Unknown orbital label '$c'")
    return i-1
end

@inline function label_from_l(l::Integer)::Char
    0 ≤ l ≤ length(L_QUANTUM_LABELS)-1 || error("Invalid l=$l")
    return L_QUANTUM_LABELS[l+1]
end

@inline function sigma_from_suffix(s::AbstractString)::Int
    su = uppercase(strip(s))
    su == "UP"   && return 1
    su == "DOWN" && return 2
    su == "↑"    && return 1
    su == "↓"    && return 2
    error("Unknown spin label '$s'")
end

@inline function suffix_from_sigma(σ::Int; style=:arrow)
    σ ∈ (1,2) || error("σ must be 1 or 2")
    if style == :arrow
        return SPIN_LABELS[σ]
    elseif style == :word
        return σ==1 ? "UP" : "DOWN"
    else
        error("style must be :arrow or :word")
    end
end

# -------------------------
# string -> tuple
# -------------------------
"""
    parse_shell(s::AbstractString)

Parse an electronic shell label.

Accepted formats
----------------
- "2s", "10d"
- "1sUP", "1sDOWN"
- "1s↑", "1s↓"

Returns
-------
- (n, l)       if no spin is specified
- (n, l, σ)    if spin is specified (σ=1 → ↑ / UP, σ=2 → ↓ / DOWN)

Notes
-----
- Uses the convention l = 0,1,2,... ↔ s,p,d,...
- Enforces n > l

Examples
--------
parse_shell("2s")      == (2,0)
parse_shell("1sDOWN")  == (1,0,2)
parse_shell("1s↓")     == (1,0,2)
"""
function parse_shell(s::AbstractString)
    t = strip(s)
    m = match(r"^(\d+)\s*([spdfghi])\s*(UP|DOWN|↑|↓)?$"i, t)
    m === nothing && error("Cannot parse '$s'")

    n = parse(Int, m.captures[1])
    l = l_from_label(lowercase(m.captures[2])[1])

    n > l || error("Convention violated: need n>l, got n=$n l=$l")

    if m.captures[3] === nothing
        return (n,l)
    else
        σ = sigma_from_suffix(m.captures[3])
        return (n,l,σ)
    end
end

# -------------------------
# tuple -> string
# -------------------------
"""
    shell_string(t::Tuple; style=:arrow)

Convert quantum numbers to an electronic shell label.

Arguments
---------
- t : (n,l) or (n,l,σ)
- style : :arrow (default) → "↑","↓"
          :word  → "UP","DOWN"

Returns
-------
String representation:
- (n,l)     → "2s"
- (n,l,σ)   → "2s↑" or "2sDOWN"

Notes
-----
- Uses the convention l = 0,1,2,... ↔ s,p,d,...
- Enforces n > l

Examples
--------
shell_string((2,0))                == "2s"
shell_string((1,0,2))              == "1s↓"
shell_string((1,0,2);style=:word)  == "1sDOWN"
"""
function shell_string(t::Tuple; style=:arrow)
    if length(t)==2
        n,l = t
        n>l || error("Need n>l")
        return string(n, label_from_l(l))
    elseif length(t)==3
        n,l,σ = t
        n>l || error("Need n>l")
        return string(n, label_from_l(l),
                      suffix_from_sigma(σ; style))
    else
        error("Tuple must be (n,l) or (n,l,σ)")
    end
end

# -------------------------
# Examples
# -------------------------
# parse_shell("2s")     -> (2,0)
# parse_shell("1sUP")   -> (1,0,1)
# parse_shell("1s↓")    -> (1,0,2)
# shell_string((2,0))           -> "2s"
# shell_string((1,0,2))         -> "1s↓"
# shell_string((1,0,2);style=:word) -> "1sDOWN"
