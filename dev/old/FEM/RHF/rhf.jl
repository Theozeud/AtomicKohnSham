using KohnShamResolution
using LinearAlgebra
using LambertW
using Plots

# Solve the eigen problem
function _eigen(basis, ρ, Rmin, Rmax; l = 0, ϵ = 1, z = 1)
    deriv_basis = deriv(basis)
    @show A = Symmetric(mass_matrix(deriv_basis))
    @show M₀ = Symmetric(mass_matrix(basis))
    @show M₋₁ = Symmetric(weight_mass_matrix(basis, -1))
    # Compute Hartree matrix
    @show F = weight_mass_vector(basis, x->4π * ρ(x) * x)
    @show Cᵨ = 4π * integrate(ρ * Monomial(2), Rmin, Rmax)
    @show Hart = vectorweight_mass_matrix(basis, A\F, Monomial(-1))+Cᵨ/(Rmax-Rmin) * M₀
    if l == 0
        H = T(0.5) * A - z .* M₋₁ + ϵ * Hart
        return eigen(H, M₀)
    else
        @show M₋₂ = Symmetric(weight_mass_matrix(basis, -2))
        @show H = T(0.5) * A - z .* M₋₁ + T(0.5) * l * (l+1) * M₋₂ + ϵ * Hart
        return eigen(H, M₀)
    end
end

function _build_density(basis, vects, z = 1)
    vect1 = vects[:, 1, 1] #1s
    vect2 = vects[:, 1, 2] #2s
    vect3 = vects[:, 2, 2] #2p
    tmp1 = build_on_basis(basis, vect1)
    tmp1 = tmp1 / sqrt(scalar_product(tmp1, tmp1, Rmin, Rmax))
    tmp1 = 1/4π * tmp1 * tmp1 * Monomial(-2) * z * 2
    tmp2 = build_on_basis(basis, vect2)
    tmp2 = tmp2 / sqrt(scalar_product(tmp2, tmp2, Rmin, Rmax))
    tmp2 = 1/4π * tmp2 * tmp2 * Monomial(-2) * z * 2
    tmp3 = build_on_basis(basis, vect3)
    tmp3 = tmp3 / sqrt(scalar_product(tmp3, tmp3, Rmin, Rmax))
    tmp3 = 1/4π * tmp3 * tmp3 * Monomial(-2) * z
    return tmp1 + tmp2 + tmp3
end

function _relax(ρ, ρprev, t)
    ρ*t + (1-t)*ρprev
end

function _iter(basis, Rmin, Rmax; l = 0, ϵ = 1, tol = 1e-2, maxiter = 10, oda = 0.8, z = 1)
    ρprev = Monomial(0, 0)
    crit = 10*tol
    iter = 0
    arraycrit = []
    arrayvals = []
    vals = zeros(length(basis), l+1)
    vects = zeros(length(basis), length(basis), l+1)
    while crit > tol && iter < maxiter
        for k in 0:l
            @time "k = $k, crit = $crit" vals[
                :, k + 1],
            vects[:, :, k + 1] = _eigen(basis, ρprev, Rmin, Rmax; l = k, ϵ = ϵ, z = z)
        end
        println("build ddensity")
        ρ = _relax(_build_density(basis, vects, z), ρprev, oda)
        println("Compute criteria")
        crit = _crit(ρ, ρprev, Rmin, Rmax)
        push!(arraycrit, crit)
        push!(arrayvals, min(vals...))
        ρprev = ρ
        println("nextiter")
        iter += 1
    end
    ρprev, arrayvals, arraycrit
end

function _crit(ρ, ρprev, Rmin, Rmax)
    X = range(Rmin, Rmax, 100)
    sqrt((Rmax-Rmin)*sum(abs.(ρ.(X[2:end]) .- ρprev.(X[2:end])) .^ 2)/(length(X)-1))
end

T = Float64
Rmin = T(0)
z = 5
@show Rmax = Int(round(-lambertw(-1e-16/sqrt(z), -1)/z))
Nmesh = 60
m = linmesh(Rmin, Rmax, Nmesh; T = T)
basis = ShortP1IntLegendreBasis(
    m, T; left = false, right = false, normalize = true, ordermin = 2, ordermax = 2)
basis_p1 = ShortP1Basis(m, T; left = false, right = false, normalize = true)
-z^2/(2)

ρ, arrayvals,
arraycrit = _iter(
    basis, Rmin, Rmax; l = 1, ϵ = 1, tol = 1e-4, maxiter = 3, oda = 0.8, z = z)
