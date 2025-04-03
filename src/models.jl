# Modelisation of exchange correlation

abstract type ExchangeCorrelation end

struct NoExchangeCorrelation <: ExchangeCorrelation end

exc(::NoExchangeCorrelation, rho) = zero(rho)
vxc(::NoExchangeCorrelation, rho) = zero(rho)

isthereExchangeCorrelation(::ExchangeCorrelation) = true
isthereExchangeCorrelation(::NoExchangeCorrelation) = false

isLSDA(::ExchangeCorrelation) = false

# SlaterXα MODEL

struct SlaterXα <: ExchangeCorrelation end

exc(::SlaterXα, rho::Real) = - 3/4 * (3/π)^(1/3) * rho^(4/3)
vxc(::SlaterXα, rho::Real) = - (3/π)^(1/3) * rho^(1/3)

exc(xa::SlaterXα, rhoUP::Real, rhoDOWN::Real) = 0.5 * (exc(xa,rhoUP) + exc(xa,rhoDOWN))
vxcUP(xa::SlaterXα, rhoUP::Real, rhoDOWN::Real) = 0.5 *vxc(xa,rhoUP)
vxcDOWN(xa::SlaterXα, rhoUP::Real, rhoDOWN::Real)  = 0.5 *vxc(xa,rhoDOWN)

## LSDA MODEL

struct LSDA <: ExchangeCorrelation end

isLSDA(::LSDA) = true

exc(lsda::LSDA, rhoUP, rhoDOWN) = ex(lsda, rhoUP,rhoDOWN) + ec(lsda, rhoUP,rhoDOWN)
vxcUP(lsda::LSDA, rhoUP, rhoDOWN) = vxUP(lsda, rhoUP) + vcUP(lsda, rhoUP, rhoDOWN)
vxcDOWN(lsda::LSDA, rhoUP, rhoDOWN) = vxDOWN(lsda, rhoDOWN) + vcDOWN(lsda, rhoUP, rhoDOWN)

ex(::LSDA, rhoUP,rhoDOWN) =  -3/4 * (6/π)^(1/3) * (rhoUP^(4/3) + rhoDOWN^(4/3))
vxUP(::LSDA, rhoUP) = -(6/π * rhoUP )^(1/3)
vxDOWN(::LSDA, rhoDOWN) = -(6/π * rhoDOWN )^(1/3)

@inline _rho(rhoUP, rhoDOWN) = rhoDOWN + rhoUP
@inline relative_spin_polarization(rhoUP, rhoDOWN, rho) = (rhoUP - rhoDOWN)/rho
@inline density_parameter(rho) = (3/(4π*rho))^(1/3)


function G(rs, p, A, α₁, β₁, β₂, β₃, β₄)
    # Spin interpolation formula
    tmp = 2*A*(β₁*rs^(1/2) + β₂*rs + β₃*rs^(3/2) + β₄*rs^(1+p))
    -2*A*(1+α₁*rs) * (log(tmp +1) - log(tmp))
end

function εcPW(rs, ξ)
    if abs(ξ) < eps(ξ)
        εcPW0 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)     # Correlation energy densioty for ξ = 0
        return εcPW0
    else
        εcPW0 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)     # Correlation energy densioty for ξ = 0
        εcPW1 =  G(rs, LSDA_CORRELATION_PARAMETERS[:,2]...)     # Correlation energy densioty for ξ = 1
        αc    = -G(rs, LSDA_CORRELATION_PARAMETERS[:,3]...)     # Spin stiffness
        fξ    = f(ξ)                                            # f(ξ)
        ξ4    = ξ^4                                             # ξ^4
        return εcPW0 + αc * (1 - ξ4) * fξ * invd2f0 + (εcPW1 - εcPW0) * ξ4 * fξ
    end
end

function ec(::LSDA, rhoUP, rhoDOWN)
    rho = _rho(rhoUP, rhoDOWN)                                 # Density
    ξ   = relative_spin_polarization(rhoDOWN, rhoUP, rho)      # Relative spin polarization
    rs  = density_parameter(rho)                               # Density parameter
    return rho * εcPW(rs, ξ)
end


@fastmath function vcUP(::LSDA, rhoUP, rhoDOWN)
    rho = _rho(rhoUP, rhoDOWN)                                  # Density
    ξ   = relative_spin_polarization(rhoUP, rhoDOWN, rho)       # Relative spin polarization
    rs  = density_parameter(rho)                                # Density parameter
    fξ    = f(ξ)                                            
    ξ4    = ξ^4

    εcPW0, ∂rs∂εcPW0 = G∂G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)     # Correlation energy densioty for ξ = 0
    εcPW1,∂rs∂εcPW1  = G∂G(rs, LSDA_CORRELATION_PARAMETERS[:,2]...)     # Correlation energy densioty for ξ = 1
    mαc,mdαc         = G∂G(rs, LSDA_CORRELATION_PARAMETERS[:,3]...)     # Minus Spin stiffness

    ΔεcPW10  = εcPW1 - εcPW0

    c1 = (1 - ξ4) * invd2f0
    c2 = ξ4 * fξ
    c3 = f(ξ) * c1

    εcPWξ  =  εcPW0 - mαc * c3 + ΔεcPW10 * c2
                                                  
    ∂rs∂εcPW =  ∂rs∂εcPW0 + (∂rs∂εcPW1 - ∂rs∂εcPW0) * c2 - mdαc * c3
    
    ∂ξ∂εcPW  = 4*ξ^3 *fξ * (ΔεcPW10 + mαc * invd2f0) + derivf(ξ) * (ξ4 * ΔεcPW10 - mαc * c1)

    return  εcPWξ - (rs/3 * ∂rs∂εcPW + (ξ - 1) * ∂ξ∂εcPW)
end

@fastmath function vcDOWN(::LSDA, rhoUP, rhoDOWN)
    rho = _rho(rhoUP, rhoDOWN)                                  # Density
    ξ   = relative_spin_polarization(rhoUP, rhoDOWN, rho)       # Relative spin polarization
    rs  = density_parameter(rho)                                # Density parameter
    fξ    = f(ξ)                                            
    ξ4    = ξ^4

    εcPW0, ∂rs∂εcPW0 = G∂G(rs, LSDA_CORRELATION_PARAMETERS[:,1]...)     # Correlation energy densioty for ξ = 0
    εcPW1,∂rs∂εcPW1  = G∂G(rs, LSDA_CORRELATION_PARAMETERS[:,2]...)     # Correlation energy densioty for ξ = 1
    mαc,mdαc         = G∂G(rs, LSDA_CORRELATION_PARAMETERS[:,3]...)     # Minus Spin stiffness

    ΔεcPW10  = εcPW1 - εcPW0

    c1 = (1 - ξ4) * invd2f0
    c2 = ξ4 * fξ
    c3 = f(ξ) * c1

    εcPWξ  =  εcPW0 - mαc * c3 + ΔεcPW10 * c2
                                                  
    ∂rs∂εcPW =  ∂rs∂εcPW0 + (∂rs∂εcPW1 - ∂rs∂εcPW0) * c2 - mdαc * c3
    
    ∂ξ∂εcPW  = 4*ξ^3 *fξ * (ΔεcPW10 + mαc * invd2f0) + derivf(ξ) * (ξ4 * ΔεcPW10 - mαc * c1)

    return εcPWξ - (rs/3 * ∂rs∂εcPW + (ξ + 1) * ∂ξ∂εcPW)
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

@fastmath function G∂G(rs, p, A, α₁, β₁, β₂, β₃, β₄)
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
    tmp     = (log(Q1 +1) - log(Q1))
    # Final quantity
    G  =  Q0 * tmp
    ∂G =  -2 * A * α₁ * tmp - (Q0*derivQ1)/(Q1^2 + Q1)
    return -G,∂G
end
 

f(ξ) = ((1 + ξ)^(4/3) + (1 - ξ)^(4/3) - 2)/(2^(4/3) - 2)
derivf(ξ) = 4/3 * ((1 + ξ)^(1/3) - (1 - ξ)^(1/3))/(2^(4/3) - 2)

const invd2f0  = (9 * (2^(1/3) -1))/4       # 1/f''(0)



# KohnSham Model

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
isLSDA(km::KohnShamExtended) = isLSDA(km.exc)

ReducedHartreeFock(z::Real, N::Real) = KohnShamExtended(z = z, N = N)
SlaterXα(z::Real, N::Real) = KohnShamExtended(z = z, N = N, exc = SlaterXα())
LSDA(z::Real, N::Real) = KohnShamExtended(z = z, N = N, exc = LSDA())