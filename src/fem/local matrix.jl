#--------------------------------------------------------------------
#                           FILL LOCAL MATRIX
#--------------------------------------------------------------------

function _fill_local_matrix!(K::AbstractArray,
                            method::FEMIntegrationMethod,
                            weight::AbstractWeight,
                            eldata::ElementData,
                            ps::PolySet,
                            basis::FEMBasis)

    @unpack a, b = eldata

    !has_singularity(weight, (a, b)) ||
                return fill_local_matrix_withsingularity!(K, method, weight, eldata, ps, basis)
    return fill_local_matrix!(K, method, weight, eldata, ps, basis)
end

#--------------------------------------------------------------------
#                           SINGULARITY INTEGRATION
#--------------------------------------------------------------------
#=
    --> Check singularity
    si oui -> on transforme l'intégration
     ∫ₐᵇ 1/x (Qi Qj)∘ϕ(x) dx
     singularité = 0 dans [a,b]
    peut on la résoudre :
    oui si Qi∘ϕ(0)= 0
    Dans ce cas, on a ϕ(x) = dϕ*x + ϕ(0)
    Or on peut écrie Qi(x) = (x - ϕ(0)) * Pi(x)
    et donc Qi∘ϕ(x)/x = (dϕ*x + ϕ(0) - ϕ(0))/x * P(ϕ(x)) = dϕ * Pi(ϕ(x))
    ∫ₐᵇ 1/x (Qi Qj)∘ϕ(x) dx  = dϕ * ∫ₐᵇ (Pi Qj)∘ϕ(x) dx
    Pour trouver Pi, on doit résoudre un problème linéaire cQ = A * cP
    où A est une matrice de taille deg(Q)*deg(Q), donc temps négligeable
=#


function fill_local_matrix_withsingularity!(K::AbstractArray,
                                            method::FEMIntegrationMethod,
                                            ::InvX,
                                            eldata::ElementData,
                                            ps::PolySet,
                                            basis::FEMBasis)
    # FACTORIZED POLYNOMIALS
    generators_factorized = factorize(basis.generators.polynomials, eldata.ϕ[2])
    # COMPUTE PRODUCT
    if eldata.s == :M
        newps = mul(generators_factorized, basis.generators.polynomials)
        fill_local_matrix!(K, method, NoWeight(), eldata, newps, basis)
        K .*= eldata.ϕ[1]
    elseif eldata.s == :T
        newps = mul(generators_factorized, basis.cache.prodMG)
        fill_local_matrix!(K, method, NoWeight(), eldata, newps, basis)
        K .*= eldata.ϕ[1]
    end
end

function fill_local_matrix_withsingularity!(K::AbstractArray,
                                            method::FEMIntegrationMethod,
                                            ::InvX2,
                                            eldata::ElementData,
                                            ps::PolySet,
                                            basis::FEMBasis)
    # FACTORIZED POLYNOMIALS
    generators_factorized = factorize(basis.generators.polynomials, eldata.ϕ[2])
    # COMPUTE PRODUCT
    if eldata.s == :M
        newps = mul(generators_factorized, generators_factorized)
        fill_local_matrix!(K, method, NoWeight(), eldata, newps, basis)
        K .*= eldata.ϕ[1]^2
    elseif eldata.s == :T
        _newps = mul(generators_factorized, generators_factorized)
        newps = mul(_newps, basis.generators.polynomials)
        fill_local_matrix!(K, method, NoWeight(), eldata, newps, basis)
        K .*= eldata.ϕ[1]^2
    end
end

#--------------------------------------------------------------------
#                           EXACT INTEGRATION
#--------------------------------------------------------------------

function fill_local_matrix!(K::AbstractArray,
                            ::ExactIntegration,
                            ::NoWeight,
                            eldata::ElementData,
                            ps::PolySet,
                            basis::FEMBasis)
    @unpack binf, bsup = eldata
    Kreshape = reshape(K,length(K),1)
    integrate!(Kreshape, ps, binf, bsup)
    Kreshape .*= eldata.invϕ[1]
end

function fill_local_matrix!(K::AbstractArray,
                            ::ExactIntegration,
                            ::InvX,
                            eldata::ElementData,
                            ps::PolySet,
                            basis::FEMBasis)
    cach = getcache(basis.cache, eldata.s)
    @unpack binf, bsup, invϕ = eldata
    @inbounds for i ∈ eachindex(cach)
        cach[i] = _integration_monome_over_deg1(i-1, invϕ[1], invϕ[2], binf, bsup)
    end
    Kreshape = reshape(K,length(K),1)
    mul!(Kreshape, ps.coeffs, cach)
    Kreshape .*= eldata.invϕ[1]
end

function fill_local_matrix!(K::AbstractArray,
                            ::ExactIntegration,
                            ::InvX2,
                            eldata::ElementData,
                            ps::PolySet,
                            basis::FEMBasis)
    cach = getcache(basis.cache, eldata.s)
    @unpack binf, bsup, invϕ = eldata
    @inbounds for i ∈ eachindex(cach)
        cach[i] = _integration_monome_over_deg2(i-1, invϕ[1], invϕ[2], binf, bsup)
    end
    Kreshape = reshape(K,length(K),1)
    mul!(Kreshape, ps.coeffs, cach)
    Kreshape .*= eldata.invϕ[1]
end


#--------------------------------------------------------------------
#                   QUADRATURE INTEGRATION
#--------------------------------------------------------------------

function fill_local_matrix!(K::AbstractArray,
                            quadra::GaussLegendre,
                            weight::FunWeight,
                            eldata::ElementData,
                            ::PolySet,
                            basis::FEMBasis)
    @unpack x, w, shiftx, fx, Qgenx = quadra
    @unpack invϕ = eldata
    shiftx .= invϕ[1].*x .+ invϕ[2]
    weight(fx, shiftx)
    @. fx  *= w
    @tensor K[i,j] = Qgenx[i,j,k] * fx[k]
    K .*= eldata.invϕ[1]
end
