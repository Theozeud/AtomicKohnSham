"""
Compute U1 = Tr(B_ℓ=1^2) bound for number of p-bound states (ℓ=1)
from a radial density rho(r) sampled on a linear grid X.

Inputs:
  X    :: AbstractVector{<:Real}   radial grid (must be sorted increasing, >= 0)
  rhoX :: AbstractVector{<:Real}   density values rho(X[i]) (>= 0 ideally)
  Z    :: Real                     nuclear charge
Keyword:
  use_positive_part :: Bool = true   if true, W = max(Z/r - V_H, 0)
                                     if false, W = Z/r - V_H (assumes V<0 everywhere)

Returns:
  U1 :: Float64
"""
function compute_U1_from_density(X, rhoX, Z; use_positive_part::Bool=true)
    @assert length(X) == length(rhoX) "X and rhoX must have same length"
    n = length(X)
    @assert n >= 3 "Need at least 3 grid points"
    @assert issorted(X) "X must be sorted increasing"
    @assert X[1] >= 0 "Grid must start at r>=0"

    # Convert to Float64 arrays (safe + fast)
    r   = Float64.(X)
    rho = Float64.(rhoX)
    Z   = Float64(Z)

    # Basic trapezoid cumulative integrator:
    # cumtrapz(f)_i ≈ ∫_{r1}^{ri} f(s) ds (with cum[1]=0)
    function cumtrapz(r, f)
        n = length(r)
        out = zeros(Float64, n)
        for i in 2:n
            dr = r[i] - r[i-1]
            out[i] = out[i-1] + 0.5 * dr * (f[i] + f[i-1])
        end
        return out
    end

    # --------- Build Q(r) and T(r) ----------
    # Q(r) = 4π ∫_0^r s^2 rho(s) ds
    fQ = 4π .* (r.^2) .* rho
    Q  = cumtrapz(r, fQ)  # Q[1]=0

    # T(r) = 4π ∫_r^∞ s rho(s) ds  (approx ∞ by last grid point)
    fT = 4π .* r .* rho
    # Reverse cumulative trapezoid: T[i] ≈ ∫_{r[i]}^{r[end]} fT(s) ds
    T = zeros(Float64, n)
    for i in (n-1):-1:1
        dr = r[i+1] - r[i]
        T[i] = T[i+1] + 0.5 * dr * (fT[i] + fT[i+1])
    end
    T[n] = 0.0

    # --------- Compute Hartree potential V_H and W ----------
    # V_H(r) = Q(r)/r + T(r). Need care at r=0.
    VH = similar(r)
    W  = similar(r)

    # Handle r=0 safely: for small r, Q(r)/r ~ O(r^2), so set Q/r = 0 at r=0.
    # More accurate: if r[1]==0, we can set VH[1] = T[1], and W[1] = +Inf formally;
    # but U1 integrals include r^-2 and will be evaluated by trapezoid, so we should
    # avoid including r=0 in integrals. We'll later start integration from the first
    # strictly positive grid point.
    for i in 1:n
        if r[i] == 0.0
            VH[i] = T[i]                  # since Q(0)=0
            W[i]  = Inf                   # Z/r diverges
        else
            VH[i] = Q[i] / r[i] + T[i]
            Wi = Z / r[i] - VH[i]
            W[i] = use_positive_part ? max(Wi, 0.0) : Wi
        end
    end

    # --------- Compute U1 with stable 1D cumulative form ----------
    # U1 = (2/9) ∫ W(r) r^-2 * A(r) dr, where A(r) = ∫_0^r W(s) s^4 ds.
    fA = W .* (r.^4)

    # Find first index with r>0 to avoid r=0 singularities in r^-2
    i0 = findfirst(>(0.0), r)
    @assert i0 !== nothing "Need at least one positive r grid point"

    # Build A(r) by cumulative trapezoid starting from i0
    A = zeros(Float64, n)
    # Set A[i0] to 0 (integral from 0 to r[i0] approximated as 0 if r[i0] small)
    for i in (i0+1):n
        dr = r[i] - r[i-1]
        A[i] = A[i-1] + 0.5 * dr * (fA[i] + fA[i-1])
    end

    # Integrand for outer integral
    integrand = zeros(Float64, n)
    for i in i0:n
        integrand[i] = W[i] * A[i] / (r[i]^2)
    end

    # Outer trapezoid integral over [r[i0], r[end]]
    I = 0.0
    for i in (i0+1):n
        dr = r[i] - r[i-1]
        I += 0.5 * dr * (integrand[i] + integrand[i-1])
    end

    U1 = (2.0 / 9.0) * I
    return U1
end

U1 = compute_U1_from_density(X, ρX, 10.0)  # Lithium: Z=3
println("U1 = ", U1)
println("Criterion: if U1 < 1 then no p bound state (n1=0).")

# trapezoidal integral on 1D grid
function trapz(x::AbstractVector{<:Real}, f::AbstractVector{<:Real})
    @assert length(x) == length(f)
    s = 0.0
    for i in 2:length(x)
        dx = float(x[i] - x[i-1])
        s += 0.5 * dx * (float(f[i]) + float(f[i-1]))
    end
    return s
end

"""
Compute I = ∫ r^3 W(r)^2 dr from rho(r) on grid X>0.
Here W(r) = V_H(rho)(r) - Z/r, with V_H the spherical Hartree potential.

Inputs:
  X    : radial grid (strictly positive, increasing)
  rhoX : rho(X[i])
  Z    : nuclear charge (for neutral atom, should satisfy 4π∫ r^2 rho = Z)

Returns:
  I :: Float64
Also returns W vector if you want it.
"""
function I_r3W2_via_W(X, rhoX, Z)
    @assert length(X) == length(rhoX)
    @assert issorted(X)
    @assert X[1] > 0

    r   = Float64.(X)
    rho = Float64.(rhoX)
    Z   = Float64(Z)
    n   = length(r)

    # Q(r) = 4π ∫_0^r s^2 rho(s) ds (approx integral from r[1] to r[i])
    fQ = 4π .* (r.^2) .* rho
    Q = zeros(Float64, n)
    for i in 2:n
        dr = r[i] - r[i-1]
        Q[i] = Q[i-1] + 0.5 * dr * (fQ[i] + fQ[i-1])
    end
    # If your grid does not start extremely close to 0, you may want to add a tiny
    # correction for [0, r[1]]; often negligible if r[1] is small.

    # T(r) = 4π ∫_r^∞ s rho(s) ds (approx to r[end])
    fT = 4π .* r .* rho
    T = zeros(Float64, n)
    for i in (n-1):-1:1
        dr = r[i+1] - r[i]
        T[i] = T[i+1] + 0.5 * dr * (fT[i] + fT[i+1])
    end

    # V_H(r) = Q(r)/r + T(r)
    VH = Q ./ r .+ T

    # W(r) = V_H(r) - Z/r  (this is ≤0 in your screenshot convention)
    W = VH .- (Z ./ r)

    # I = ∫ r^3 W(r)^2 dr
    integrand = (r.^3) .* (W.^2)
    I = trapz(r, integrand)
    return I, W
end

I, W = I_r3W2_via_W(X, ρX, 10)
println("I = ∫ r^3 W^2 = ", 2/15*I)
