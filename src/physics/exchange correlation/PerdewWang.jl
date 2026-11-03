struct PerdewWang <: BuiltinFunctional
    identifier::Symbol
    n_spin::Int
    name::String
    family::Symbol
    function PerdewWang(n_spin::Int)
        new(:lda_c_pw, n_spin, "Perdew-Wang Correlation", :lda)
    end
end

# Perdew & Wang, Phys. Rev. B 45, 13244 (1992), LDA correlation.
# G(rs; A,α1,β1,β2,β3,β4) is the shared rational/log form used for εc(rs,0),
# εc(rs,1) and the spin-stiffness -αc(rs); parameter sets P0/P1/Pα below are
# Table I of that paper.
function G(rs::Real, A, α1, β1, β2, β3, β4)
    Q0 = -2A * (1 + α1 * rs)
    Q1 = 2A * (β1 * sqrt(rs) + β2 * rs + β3 * rs^(3//2) + β4 * rs^2)
    Q0 * log(1 + 1 / Q1)
end

function dG_drs(rs::Real, A, α1, β1, β2, β3, β4)
    Q0 = -2A * (1 + α1 * rs)
    dQ0 = -2A * α1
    Q1 = 2A * (β1 * sqrt(rs) + β2 * rs + β3 * rs^(3//2) + β4 * rs^2)
    dQ1 = 2A * (β1 / (2sqrt(rs)) + β2 + 3//2 * β3 * sqrt(rs) + 2β4 * rs)
    dQ0 * log(1 + 1 / Q1) - Q0 * dQ1 / (Q1 * (Q1 + 1))
end

const PW92_P0 = (0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294)
const PW92_P1 = (0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517)
const PW92_Pα = (0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)

rs_of_ρ(ρ::Real) = (3 / (4π * ρ))^(1//3)

eval_zk(::PerdewWang, ρ::Real) = iszero(ρ) ? zero(ρ) : G(rs_of_ρ(ρ), PW92_P0...)

function eval_vrho(::PerdewWang, ρ::Real)
    iszero(ρ) && return zero(ρ)
    rs = rs_of_ρ(ρ)
    ec = G(rs, PW92_P0...)
    decdrs = dG_drs(rs, PW92_P0...)
    ec - rs / 3 * decdrs
end

# Spin-polarized PW92: the ζ-interpolation between the ζ=0 and ζ=1 correlation
# energies, weighted by the spin-stiffness αc(rs) = -G(rs, Pα...) through the
# spin-scaling function f(ζ) = [(1+ζ)^(4/3)+(1-ζ)^(4/3)-2]/(2^(4/3)-2), whose
# curvature at the origin is f''(0) = 4/(9(2^(1/3)-1)) (PW92_F2P0 below).
#
#   εc(rs,ζ) = εc(rs,0) + αc(rs)[f(ζ)/f''(0)](1-ζ⁴) + [εc(rs,1)-εc(rs,0)]f(ζ)ζ⁴
#
# Validated against Libxc's lda_c_pw (n_spin=2) to ~1e-9 relative error across
# ζ ∈ [0,1], including the ζ=0 and ζ=1 edge cases (machine precision there).
const PW92_F2P0 = 4 / (9 * (2^(1//3) - 1))
pw92_fζ(ζ::Real) = ((1 + ζ)^(4//3) + (1 - ζ)^(4//3) - 2) / (2^(4//3) - 2)
pw92_dfζ(ζ::Real) = 4//3 * ((1 + ζ)^(1//3) - (1 - ζ)^(1//3)) / (2^(4//3) - 2)

function pw92_ec(rs::Real, ζ::Real)
    ec0 = G(rs, PW92_P0...)
    ec1 = G(rs, PW92_P1...)
    αc = -G(rs, PW92_Pα...)
    f = pw92_fζ(ζ)
    ec0 + αc * (f / PW92_F2P0) * (1 - ζ^4) + (ec1 - ec0) * f * ζ^4
end

function pw92_decdrs(rs::Real, ζ::Real)
    dec0 = dG_drs(rs, PW92_P0...)
    dec1 = dG_drs(rs, PW92_P1...)
    dαc = -dG_drs(rs, PW92_Pα...)
    f = pw92_fζ(ζ)
    dec0 + dαc * (f / PW92_F2P0) * (1 - ζ^4) + (dec1 - dec0) * f * ζ^4
end

function pw92_decdζ(rs::Real, ζ::Real)
    ec0 = G(rs, PW92_P0...)
    ec1 = G(rs, PW92_P1...)
    αc = -G(rs, PW92_Pα...)
    f = pw92_fζ(ζ)
    df = pw92_dfζ(ζ)
    αc / PW92_F2P0 * (df * (1 - ζ^4) - 4ζ^3 * f) + (ec1 - ec0) * (df * ζ^4 + 4ζ^3 * f)
end

function eval_zk(::PerdewWang, ρup::Real, ρdown::Real)
    ρ = ρup + ρdown
    iszero(ρ) && return zero(ρ)
    ζ = (ρup - ρdown) / ρ
    pw92_ec(rs_of_ρ(ρ), ζ)
end

function _pw92_vup_vdown(ρup::Real, ρdown::Real)
    ρ = ρup + ρdown
    ζ = (ρup - ρdown) / ρ
    rs = rs_of_ρ(ρ)
    ec = pw92_ec(rs, ζ)
    decdrs = pw92_decdrs(rs, ζ)
    decdζ = pw92_decdζ(rs, ζ)
    common = ec - rs / 3 * decdrs
    common + (1 - ζ) * decdζ, common - (1 + ζ) * decdζ
end
eval_vrho_up(func::PerdewWang, ρup::Real, ρdown::Real) = iszero(ρup + ρdown) ? zero(ρup) : _pw92_vup_vdown(ρup, ρdown)[1]
eval_vrho_down(func::PerdewWang, ρup::Real, ρdown::Real) = iszero(ρup + ρdown) ? zero(ρdown) : _pw92_vup_vdown(ρup, ρdown)[2]

# ---------------------------------------------------------------------------
# Second derivatives (exchange-correlation kernel), Libxc `v2rho2` convention:
# d²(ρεc)/dρ² unpolarized; (↑↑, ↑↓, ↓↓) spin-polarized.
# ---------------------------------------------------------------------------
function d2G_drs2(rs::Real, A, α1, β1, β2, β3, β4)
    Q0 = -2A * (1 + α1 * rs)
    dQ0 = -2A * α1
    Q1 = 2A * (β1 * sqrt(rs) + β2 * rs + β3 * rs^(3//2) + β4 * rs^2)
    dQ1 = 2A * (β1 / (2sqrt(rs)) + β2 + 3//2 * β3 * sqrt(rs) + 2β4 * rs)
    d2Q1 = 2A * (-β1 / (4rs^(3//2)) + 3//4 * β3 / sqrt(rs) + 2β4)
    P = Q1 * (Q1 + 1)
    dL = -dQ1 / P
    d2L = -d2Q1 / P + dQ1^2 * (2Q1 + 1) / P^2
    2dQ0 * dL + Q0 * d2L
end

function eval_v2rho2(::PerdewWang, ρ::Real)
    iszero(ρ) && return zero(ρ)
    rs = rs_of_ρ(ρ)
    d1 = dG_drs(rs, PW92_P0...)
    d2 = d2G_drs2(rs, PW92_P0...)
    -rs / (3ρ) * (2//3 * d1 - rs / 3 * d2)
end

pw92_d2fζ(ζ::Real) = 4//9 * ((1 + ζ)^(-2//3) + (1 - ζ)^(-2//3)) / (2^(4//3) - 2)

"""
    pw92_derivatives(rs, ζ) -> (ε, ε_rs, ε_rsrs, ε_ζ, ε_ζζ, ε_rsζ)

PW92 correlation energy per particle `ε(rs, ζ)` and its first and second partial
derivatives, from `ε = ec0 + αc g(ζ) + (ec1 - ec0) h(ζ)` with
`g = f(1 - ζ⁴)/f''(0)` and `h = f ζ⁴`.
"""
function pw92_derivatives(rs::Real, ζ::Real)
    e0, e1, a = G(rs, PW92_P0...), G(rs, PW92_P1...), -G(rs, PW92_Pα...)
    de0, de1, da = dG_drs(rs, PW92_P0...), dG_drs(rs, PW92_P1...), -dG_drs(rs, PW92_Pα...)
    d2e0, d2e1, d2a = d2G_drs2(rs, PW92_P0...), d2G_drs2(rs, PW92_P1...),
                      -d2G_drs2(rs, PW92_Pα...)
    f, df, d2f = pw92_fζ(ζ), pw92_dfζ(ζ), pw92_d2fζ(ζ)
    g = f * (1 - ζ^4) / PW92_F2P0
    dg = (df * (1 - ζ^4) - 4ζ^3 * f) / PW92_F2P0
    d2g = (d2f * (1 - ζ^4) - 8ζ^3 * df - 12ζ^2 * f) / PW92_F2P0
    h = f * ζ^4
    dh = df * ζ^4 + 4ζ^3 * f
    d2h = d2f * ζ^4 + 8ζ^3 * df + 12ζ^2 * f
    ε = e0 + a * g + (e1 - e0) * h
    ε_rs = de0 + da * g + (de1 - de0) * h
    ε_rsrs = d2e0 + d2a * g + (d2e1 - d2e0) * h
    ε_ζ = a * dg + (e1 - e0) * dh
    ε_ζζ = a * d2g + (e1 - e0) * d2h
    ε_rsζ = da * dg + (de1 - de0) * dh
    ε, ε_rs, ε_rsrs, ε_ζ, ε_ζζ, ε_rsζ
end

"""
    eval_v2rho2(::PerdewWang, ρup, ρdown) -> (f↑↑, f↑↓, f↓↓)

Spin-polarized PW92 kernel `∂v_σ/∂ρ_σ'`, by the chain rule through
`rs(ρ)` and `ζ = (ρ↑ - ρ↓)/ρ`. `ζ` is kept at a distance `100eps` from `±1`,
where `f''(ζ)` diverges, so that a fully polarized density gives large finite
values instead of `Inf` (they multiply density variations that vanish there).
"""
function eval_v2rho2(::PerdewWang, ρup::Real, ρdown::Real)
    ρ = ρup + ρdown
    iszero(ρ) && return (zero(ρ), zero(ρ), zero(ρ))
    ζmax = 1 - 100eps(typeof(ρ))
    ζ = clamp((ρup - ρdown) / ρ, -ζmax, ζmax)
    rs = rs_of_ρ(ρ)
    _, ε_rs, ε_rsrs, ε_ζ, ε_ζζ, ε_rsζ = pw92_derivatives(rs, ζ)
    dc_drs = 2//3 * ε_rs - rs / 3 * ε_rsrs
    dc_dζ = ε_ζ - rs / 3 * ε_rsζ
    dvup_drs = dc_drs + (1 - ζ) * ε_rsζ
    dvdn_drs = dc_drs - (1 + ζ) * ε_rsζ
    dvup_dζ = dc_dζ - ε_ζ + (1 - ζ) * ε_ζζ
    dvdn_dζ = dc_dζ - ε_ζ - (1 + ζ) * ε_ζζ
    drs = -rs / (3ρ)
    dζ_up, dζ_dn = (1 - ζ) / ρ, -(1 + ζ) / ρ
    fuu = dvup_drs * drs + dvup_dζ * dζ_up
    fud = dvup_drs * drs + dvup_dζ * dζ_dn
    fdd = dvdn_drs * drs + dvdn_dζ * dζ_dn
    (fuu, fud, fdd)
end
