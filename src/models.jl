#####################################################################
#                        EXCHANGE-CORRELATION
#####################################################################

abstract type ExchangeCorrelation end

struct NoExchangeCorrelation <: ExchangeCorrelation end

exc(::NoExchangeCorrelation, rho) = zero(rho)
vxc(::NoExchangeCorrelation, rho) = zero(rho)

isthereExchangeCorrelation(::ExchangeCorrelation) = true
isthereExchangeCorrelation(::NoExchangeCorrelation) = false

#####################################################################
#                            SlaterXα
#####################################################################

struct SlaterXα <: ExchangeCorrelation end

exc(::SlaterXα, rho::Real) = - 3/4 * (3/π)^(1/3) * rho^(4/3)
vxc(::SlaterXα, rho::Real) = - (3/π)^(1/3) * rho^(1/3)

exc(xa::SlaterXα, rhoUP::Real, rhoDOWN::Real) = 0.5 * (exc(xa,rhoUP) + exc(xa,rhoDOWN))
vxcUP(xa::SlaterXα, rhoUP::Real, rhoDOWN::Real) = 0.5 *vxc(xa,rhoUP)
vxcDOWN(xa::SlaterXα, rhoUP::Real, rhoDOWN::Real)  = 0.5 *vxc(xa,rhoDOWN)

#####################################################################
#                           Perdew-Wang 92
#####################################################################

struct LSDA <: ExchangeCorrelation end

exc(lsda::LSDA, rhoUP::Real, rhoDOWN::Real) = ex(lsda, rhoUP,rhoDOWN) + ec(lsda, rhoUP,rhoDOWN)
vxcUP(lsda::LSDA, rhoUP::Real, rhoDOWN::Real) = vxUP(lsda, rhoUP) + vcUP(lsda, rhoUP, rhoDOWN)
vxcDOWN(lsda::LSDA, rhoUP::Real, rhoDOWN::Real) = vxDOWN(lsda, rhoDOWN) + vcDOWN(lsda, rhoUP, rhoDOWN)

ex(::LSDA, rhoUP::Real,rhoDOWN::Real) =  -3/4 * (6/π)^(1/3) * (rhoUP^(4/3) + rhoDOWN^(4/3))
vxUP(::LSDA, rhoUP::Real) = -(6/π * rhoUP )^(1/3)
vxDOWN(::LSDA, rhoDOWN::Real) = -(6/π * rhoDOWN )^(1/3)

@inline _rho(rhoUP::Real, rhoDOWN::Real) = rhoDOWN + rhoUP
@inline relative_spin_polarization(rhoUP::Real, rhoDOWN::Real, rho) = (rhoUP - rhoDOWN)/rho
@inline density_parameter(rho::Real) = (3/(4π*rho))^(1/3)

function G(rs::Real, p::Real, A::Real, α₁::Real, β₁::Real, β₂::Real, β₃::Real, β₄::Real)
    # Spin interpolation formula
    tmp = 2*A*(β₁*rs^(1/2) + β₂*rs + β₃*rs^(3/2) + β₄*rs^(1+p))
    -2*A*(1+α₁*rs) * log(1 +1/tmp)
end

function εcPW(rs::Real, ξ::Real)
    if abs(ξ) < eps(ξ)
        @views view_params1 = LSDA_CORRELATION_PARAMETERS[:,1]
        εcPW0 =  G(rs, view_params1...)     # Correlation energy densioty for ξ = 0
        return εcPW0
    else
        @views view_params1 = LSDA_CORRELATION_PARAMETERS[:,1]
        @views view_params2 = LSDA_CORRELATION_PARAMETERS[:,2]
        @views view_params3 = LSDA_CORRELATION_PARAMETERS[:,3]
        εcPW0 =  G(rs, view_params1...)    # Correlation energy densioty for ξ = 0
        εcPW1 =  G(rs, view_params2...)     # Correlation energy densioty for ξ = 1
        αc    = -G(rs, view_params3...)     # Spin stiffness
        fξ    = f(ξ)                                            # f(ξ)
        ξ4    = ξ^4                                             # ξ^4
        return εcPW0 + αc * (1 - ξ4) * fξ * invd2f0 + (εcPW1 - εcPW0) * ξ4 * fξ
    end
end

function ec(::LSDA, rhoUP::Real, rhoDOWN::Real)
    rho = _rho(rhoUP, rhoDOWN)                                 # Density
    ξ   = relative_spin_polarization(rhoUP, rhoDOWN, rho)      # Relative spin polarization
    rs  = density_parameter(rho)                               # Density parameter
    return rho * εcPW(rs, ξ)
end


function vcUP(lsda::LSDA, rhoUP::Real, rhoDOWN::Real)
    _vc(lsda, rhoUP, rhoDOWN, 1)
end

function vcDOWN(lsda::LSDA, rhoUP::Real, rhoDOWN::Real)
    _vc(lsda, rhoUP, rhoDOWN, -1)
end

@fastmath function _vc(::LSDA, rhoUP::Real, rhoDOWN::Real, σ::Int)
    # DENSITY
    rho = _rho(rhoUP, rhoDOWN)
    
    # RELATIVE SPIN POLARIZATION
    ξ   = relative_spin_polarization(rhoUP, rhoDOWN, rho)
    
    # DENSITY PARAMETER
    rs  = density_parameter(rho)   
    
    # CONSTANTS
    fξ    = f(ξ)                                            
    ξ4    = ξ^4
    c1 = (1 - ξ4) * invd2f0
    c2 = ξ4 * fξ
    c3 = fξ * c1

    # CORRELATION ENERGY DENSITY FOR ξ = 0
    @views view_params1 = LSDA_CORRELATION_PARAMETERS[:,1]
    εcPW0, ∂εcPW0∂rs = G∂G(rs, view_params1...)

    # CORRELATION ENERGY DENSITY FOR ξ = 1
    @views view_params2 = LSDA_CORRELATION_PARAMETERS[:,2]
    εcPW1,∂εcPW1∂rs  = G∂G(rs, view_params2...)

    # MINUS SPIN STIFFNESS
    @views view_params3 = LSDA_CORRELATION_PARAMETERS[:,3]
    mαc,mdαc         = G∂G(rs, view_params3...)


    ΔεcPW10  = εcPW1 - εcPW0
    Δ∂εcPW∂rs = ∂εcPW1∂rs - ∂εcPW0∂rs
    
    εcPWξ  =  εcPW0 - mαc * c3 + ΔεcPW10 * c2
                                                 
    ∂εcPW∂rs =  ∂εcPW0∂rs + Δ∂εcPW∂rs * c2 - mdαc * c3
    
    ∂εcPW∂ξ  = 4*ξ^3 *fξ * (ΔεcPW10 + mαc * invd2f0) + derivf(ξ) * (ξ4 * ΔεcPW10 - mαc * c1)

    return  εcPWξ - rs/3 * ∂εcPW∂rs - (ξ - σ) * ∂εcPW∂ξ
end



const LSDA_CORRELATION_PARAMETERS =
      # εcPW0       # εcPW1      # -αc
    [     1.0           1.0        1.0;  # p 
     0.031091      0.015545   0.016887;  # A
      0.21370       0.20548    0.11125;  # α₁
       7.5957       14.1189     10.357;  # β₁
       3.5876        6.1977     3.6231;  # β₂
       1.6382        3.3662    0.88026;  # β₃
      0.49294       0.62517    0.49671   # β₄
    ] 

@fastmath function G∂G(rs::Real, p::Real, A::Real, α₁::Real, β₁::Real, β₂::Real, β₃::Real, β₄::Real)
    # Calculs of power of rs 
    rs12    = sqrt(rs)
    rs_12   = rs12/rs
    rs32    = rs * rs12
    rsp     = rs^p
    rsp1    = rsp * rs
    # Calculs of Intermediate quantity
    Q0      = -2*A*(1+α₁*rs)
    Q1      =  2*A*(β₁*rs12 + β₂*rs + β₃*rs32 + β₄*rsp1)
    derivQ1 = A*(β₁*rs_12 + 2*β₂ + 3*β₃*rs12 + 2*(p+1)*β₄ *rsp)
    tmp     = log(1 + 1/Q1)
    # Final quantity
    G  =  Q0 * tmp
    ∂G =  -2 * A * α₁ * tmp - (Q0*derivQ1)/(Q1^2 + Q1)
    return -G,∂G
end
 

f(ξ::Real) = ((1 + ξ)^(4/3) + (1 - ξ)^(4/3) - 2)/(2^(4/3) - 2)
derivf(ξ::Real) = 4/3 * ((1 + ξ)^(1/3) - (1 - ξ)^(1/3))/(2^(4/3) - 2)

const invd2f0  = (9 * (2^(1/3) -1))/4       # 1/f''(0)



#####################################################################
#                         MODEL STRUCTURE
#####################################################################


abstract type AbstractDFTModel end 

struct KohnShamExtended{TEXCH <: ExchangeCorrelation, TZ <: Real, TN <: Real} <: AbstractDFTModel
    z::TZ
    N::TN
    exc::TEXCH

    function KohnShamExtended(;z::Real, N::Real, exc::ExchangeCorrelation = NoExchangeCorrelation())
        new{typeof(exc), eltype(z), eltype(N)}(z,N,exc)
    end
end

isthereExchangeCorrelation(km::KohnShamExtended) = isthereExchangeCorrelation(km.exc)

ReducedHartreeFock(z::Real, N::Real) = KohnShamExtended(z = z, N = N)
SlaterXα(z::Real, N::Real) = KohnShamExtended(z = z, N = N, exc = SlaterXα())
LSDA(z::Real, N::Real) = KohnShamExtended(z = z, N = N, exc = LSDA())