#--------------------------------------------------------------------
#                         MODEL STRUCTURE
#--------------------------------------------------------------------
"""
    KSEModel(; Z::Real, N::Real, hartree::Real = 0, exc::ExchangeCorrelation = NoExchangeCorrelation())

Structure for storing the parameters of an extended Kohn–Sham model (LDA/LSDA).
Encodes physical and modeling parameters for atomic or ionic simulations.

# Arguments
- `Z::Real`: Nuclear charge (number of protons).
- `N::Real`: Number of electrons.
- `hartree::Real = 0`: Scaling factor for the Hartree term:
    - `0` disables electron-electron repulsion.
    - `1` uses the full Hartree interaction.
- `exc::ExchangeCorrelation = NoExchangeCorrelation()`: Exchange–correlation model used in the simulation.
"""
struct KSEModel{T <: Real,
    TEX,
    TCO}
    Z::T                # Charge of the nucleus
    N::T                # Number of electrons

    # Coefficient multiply to the Hartree Term :
    # 0 -> no hartree term,
    # 1-> full hartree term
    hartree::T

    exchange::TEX       # Exchange functional
    correlation::TCO    # Correlation functional

    # Spin Polarization
    # 1 -> No Polarization
    # 2 -> Polarisation
    nspin::Int

    function KSEModel(; Z::Real,
                        N::Real,
                        hartree::Real = 1,
                        ex = NoFunctional(1),
                        ec = NoFunctional(1))
        T = promote_type(typeof(Z), typeof(N), typeof(hartree))
        nspin = max(ex.n_spin, ec.n_spin)
        new{T,
            typeof(ex),
            typeof(ec)}(T(Z), T(N), T(hartree), ex, ec, nspin)
    end
end

has_exchange(model::KSEModel) = !(model.exchange isa NoFunctional)
has_correlation(model::KSEModel) = !(model.correlation isa NoFunctional)
has_exchcorr(model::KSEModel) = has_exchange(model) || has_correlation(model)

function evaluate_vrhox(model::KSEModel; rho::AbstractArray{<:Real})
    evaluate_functional(model.exchange; rho = rho, derivatives = 1)
end

function evaluate_vrhoc(model::KSEModel; rho::AbstractArray{<:Real})
    evaluate_functional(model.correlation; rho = rho, derivatives = 1)
end

function evaluate_vrhox!(
        model::KSEModel;
        rho::AbstractArray{<:Real},
        vrho::AbstractArray{<:Real})
    evaluate_functional!(model.exchange; rho = rho, vrho = vrho)
end

function evaluate_vrhoc!(
        model::KSEModel;
        rho::AbstractArray{<:Real},
        vrho::AbstractArray{<:Real})
    evaluate_functional!(model.correlation; rho = rho, vrho = vrho)
end

function evaluate_vrho!(model::KSEModel;
        rho::AbstractArray{<:Real},
        vrho::AbstractArray{<:Real},
        cache::AbstractArray{<:Real})
    evaluate_vrhoc!(model; rho = rho, vrho = cache)
    evaluate_vrhox!(model; rho = rho, vrho = vrho)
    @. vrho += cache
end

function evaluate_zk!(
        model::KSEModel;
        rho::AbstractArray{<:Real},
        zk::AbstractArray{<:Real},
        cache::AbstractArray{<:Real})
    evaluate_functional!(model.correlation; rho = rho, zk = cache)
    evaluate_functional!(model.exchange; rho = rho, zk = zk)
    @. zk += cache
end

#--------------------------------------------------------------------
#                         SHORT MODEL CONSTRUCTOR
#--------------------------------------------------------------------
"""
    RHF(; Z::Real, N::Real)

Convenience constructor for the Hartree–Fock model using the extended Kohn–Sham framework.

Creates a `KSEModel` with the specified nuclear charge `Z` and number of electrons `N`, and
disables the exchange–correlation.

# Arguments
- `Z::Real`: Nuclear charge (e.g., 2.0 for helium).
- `N::Real`: Number of electrons.

# Returns
- `KSEModel` with no exchange–correlation and no Hartree term.
"""
RHF(; Z::Real, N::Real) = KSEModel(Z = Z, N = N)

"""
    Slater(; Z::Real, N::Real)

Convenience constructor for an extended Kohn–Sham model with the Slater exchange functional.

Creates a `KSEModel` with the specified nuclear charge `Z` and number of electrons `N`,
and uses the Slater exchange-only approximation (no correlation term and no Hartree term).

# Arguments
- `Z::Real`: Nuclear charge (e.g., 10.0 for neon).
- `N::Real`: Number of electrons.

# Returns
- `KSEModel` using the Slater Xα exchange functional.
"""
function Slater(; Z::Real, N::Real, nspin::Int = 1)
    ex = BuiltinFunctional(:lda_x; nspin = nspin)
    ec = NoFunctional(nspin)
    return KSEModel(Z = Z, N = N, ex = ex, ec = ec)
end
