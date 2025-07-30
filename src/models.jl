#--------------------------------------------------------------------
#                         MODEL STRUCTURE
#--------------------------------------------------------------------
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
    TEX,
    TCO}
    z::T                # Charge of the nucleus
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
    n_spin::Int

    function KSEModel(; z::Real,
                        N::Real,
                        hartree::Real = 1,
                        ex = NoFunctional(1),
                        ec = NoFunctional(1))
        T = promote_type(typeof(z), typeof(N), typeof(hartree))
        n_spin = max(ex.n_spin, ec.n_spin)
        new{T,
            typeof(ex),
            typeof(ec)}(T(z), T(N), T(hartree), ex, ec, n_spin)
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
    RHF(; z::Real, N::Real)

Convenience constructor for the Hartree–Fock model using the extended Kohn–Sham framework.

Creates a `KSEModel` with the specified nuclear charge `z` and number of electrons `N`, and
disables the exchange–correlation.

# Arguments
- `z::Real`: Nuclear charge (e.g., 2.0 for helium).
- `N::Real`: Number of electrons.

# Returns
- `KSEModel` with no exchange–correlation and no Hartree term.
"""
RHF(; z::Real, N::Real) = KSEModel(z = z, N = N)

"""
    Slater(; z::Real, N::Real)

Convenience constructor for an extended Kohn–Sham model with the Slater exchange functional.

Creates a `KSEModel` with the specified nuclear charge `z` and number of electrons `N`,
and uses the Slater exchange-only approximation (no correlation term and no Hartree term).

# Arguments
- `z::Real`: Nuclear charge (e.g., 10.0 for neon).
- `N::Real`: Number of electrons.

# Returns
- `KSEModel` using the Slater Xα exchange functional.
"""
function Slater(; z::Real, N::Real, n_spin::Int = 1)
    ex = BuiltinFunctional(:lda_x; n_spin = n_spin)
    ec = NoFunctional(n_spin)
    return KSEModel(z = z, N = N, ex = ex, ec = ec)
end
