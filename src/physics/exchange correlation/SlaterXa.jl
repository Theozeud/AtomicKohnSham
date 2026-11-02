struct SlaterXα <: BuiltinFunctional
    identifier::Symbol
    n_spin::Int
    name::String
    family::Symbol
    function SlaterXα(n_spin::Int)
        new(:lda_x, n_spin, "Slater Exchange", :lda)
    end
end

eval_zk(::SlaterXα, ρ::Real) = - 3/4 * (3/π)^(1/3) * ρ^(1/3)
eval_vrho(::SlaterXα, ρ::Real) = - (3/π)^(1/3) * ρ^(1/3)

# Spin-polarized (LSDA) Slater exchange, via the standard spin-scaling relation
# Eₓ[ρ↑,ρ↓] = (Eₓ[2ρ↑] + Eₓ[2ρ↓])/2 applied to the spin-unpolarized functional above.
#
# eval_zk is the closed form of that relation (energy per particle):
#   εₓ(ρ↑,ρ↓) = -3/4 (6/π)^(1/3) (ρ↑^(4/3) + ρ↓^(4/3)) / (ρ↑+ρ↓)
# which reduces to eval_zk(ρ) at ρ↑=ρ↓=ρ/2, matching the unpolarized case exactly.
#
# eval_vrho_up/down are the functional derivatives δEₓ/δρ↑, δEₓ/δρ↓ of that same
# relation, which only depend on the same-spin density (no cross term for LDA
# exchange): vrho_up(ρ↑,ρ↓) = eval_vrho(2ρ↑), vrho_down(ρ↑,ρ↓) = eval_vrho(2ρ↓).
function eval_zk(func::SlaterXα, ρup::Real, ρdown::Real)
    ρ = ρup + ρdown
    iszero(ρ) && return zero(ρ)
    - 3/4 * (6/π)^(1/3) * (ρup^(4/3) + ρdown^(4/3)) / ρ
end
eval_vrho_up(func::SlaterXα, ρup::Real, ρdown::Real) = eval_vrho(func, 2ρup)
eval_vrho_down(func::SlaterXα, ρup::Real, ρdown::Real) = eval_vrho(func, 2ρdown)

# Second derivatives (exchange-correlation kernel), Libxc `v2rho2` convention:
# d²(ρεₓ)/dρ² unpolarized; (↑↑, ↑↓, ↓↓) spin-polarized. Exchange has no
# cross-spin term. The kernel diverges like ρ^(-2/3) at ρ = 0, where it is set to 0.
eval_v2rho2(::SlaterXα, ρ::Real) = iszero(ρ) ? zero(ρ) : - (3/π)^(1/3) / 3 * ρ^(-2/3)
function eval_v2rho2(func::SlaterXα, ρup::Real, ρdown::Real)
    (2eval_v2rho2(func, 2ρup), zero(ρup), 2eval_v2rho2(func, 2ρdown))
end
