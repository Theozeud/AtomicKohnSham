#--------------------------------------------------------------------
#                         MODEL STRUCTURE
#--------------------------------------------------------------------

# RAJOUTER LE TERME DE HARTREE
"""
    KSEModel(; z::Real, N::Real, hartree::Real = 0, exc::ExchangeCorrelation = NoExchangeCorrelation())

Structure for storing the parameters of an extended Kohn–Sham model (LDA/LSDA).  
Encodes physical and modeling parameters for atomic or ionic simulations.

# Arguments
- `z::Real`: Nuclear charge (number of protons).
- `N::Real`: Number of electrons.
- `hartree::Real = 0`: Scaling factor for the Hartree term:
    - `0` disables electron-electron repulsion.
    - `1` uses the full Hartree interaction.
- `exc::ExchangeCorrelation = NoExchangeCorrelation()`: Exchange–correlation model used in the simulation. 
"""
struct KSEModel{T <: Real, 
                TEXCH <: ExchangeCorrelation}
    z::T
    N::T
    hartree::T
    exc::TEXCH
    function KSEModel(; z::Real, 
                        N::Real, 
                        hartree::Real = 1, 
                        exc::ExchangeCorrelation = NoExchangeCorrelation())
        T = promote_type(z, N, Hartree)
        new{T, typeof(exc)}(T(z), T(N), T(hartree), exc)
    end
end

isthereExchangeCorrelation(km::KSEModel) = isthereExchangeCorrelation(km.exc)


#--------------------------------------------------------------------
#                         MODEL CONSTRUCTOR
#--------------------------------------------------------------------
"""
    RHF(; z::Real, N::Real)

Convenience constructor for the Hartree–Fock model using the extended Kohn–Sham framework.

Creates a `KSEModel` with the specified nuclear charge `z` and number of electrons `N`, and disables the exchange–correlation and Hartree terms, corresponding to the standard Hartree–Fock approximation.

# Arguments
- `z::Real`: Nuclear charge (e.g., 2.0 for helium).
- `N::Real`: Number of electrons.

# Returns
- `KSEModel` with no exchange–correlation and no Hartree term.
"""
RHF(; z::Real, N::Real) = KSEModel(z = z, N = N)


"""
    SlaterXα(; z::Real, N::Real)

Convenience constructor for an extended Kohn–Sham model with the Slater Xα exchange functional.

Creates a `KSEModel` with the specified nuclear charge `z` and number of electrons `N`, and uses the Slater Xα exchange-only approximation (no correlation term and no Hartree term).

# Arguments
- `z::Real`: Nuclear charge (e.g., 10.0 for neon).
- `N::Real`: Number of electrons.

# Returns
- `KSEModel` using the Slater Xα exchange functional.
"""
SlaterXα(;z::Real, N::Real) = KSEModel(z = z, N = N, exc = SlaterXα())



# MUST BE RENAMED IN SOMETHING LIKE PERDEW 
LSDA(;z::Real, N::Real) = KSEModel(z = z, N = N, exc = LSDA())