#--------------------------------------------------------------------
#                         MODEL STRUCTURE
#--------------------------------------------------------------------
"""
    KSEModel(; Z::Real, N::Real, hartree::Real = 1, ex = NoFunctional(1), ec = NoFunctional(1))

Structure for storing the parameters of an extended Kohn–Sham model (LDA/LSDA).
Encodes physical and modeling parameters for atomic or ionic simulations.

# Arguments
- `Z::Real`: Nuclear charge (number of protons).
- `N::Real`: Number of electrons.
- `hartree::Real = 1`: Scaling factor for the Hartree term:
    - `0` disables electron-electron repulsion.
    - `1` uses the full Hartree interaction.
- `ex`: Exchange functional (e.g. a Libxc `Functional`, or `NoFunctional(nspin)` to disable it).
- `ec`: Correlation functional (same convention as `ex`).
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
        # Mixed-rung XC (LDA exchange + GGA correlation or vice versa) is not
        # implemented: the GGA weak-form assembly (assemble_exc!) assumes
        # both contributions need ∇ρ or neither does. NoFunctional is exempt
        # since it contributes nothing either way.
        if !(ex isa NoFunctional) && !(ec isa NoFunctional) && is_gga(ex) != is_gga(ec)
            error("Mixed-rung exchange-correlation is not supported: exchange is " *
                  "$(is_gga(ex) ? "GGA" : "LDA") while correlation is " *
                  "$(is_gga(ec) ? "GGA" : "LDA").")
        end
        new{T,
            typeof(ex),
            typeof(ec)}(T(Z), T(N), T(hartree), ex, ec, nspin)
    end
end

has_exchange(model::KSEModel) = !(model.exchange isa NoFunctional)
has_correlation(model::KSEModel) = !(model.correlation isa NoFunctional)
has_exchcorr(model::KSEModel) = has_exchange(model) || has_correlation(model)
has_gga(model::KSEModel) = is_gga(model.exchange) || is_gga(model.correlation)

function evaluate_vrhox(model::KSEModel; rho::AbstractArray{<:Real})
    evaluate_functional(model.exchange; rho = rho, derivatives = 1)
end

function evaluate_vrhoc(model::KSEModel; rho::AbstractArray{<:Real})
    evaluate_functional(model.correlation; rho = rho, derivatives = 1)
end


# Only forward `sigma` when actually given: Libxc's own evaluate!/evaluate
# treats a `sigma` keyword's mere presence (even `nothing`) as "this argument
# was supplied" and tries `length(sigma)` on it -- so LDA call sites (sigma
# left at its `nothing` default) must not pass the keyword at all.
_sigma_kwargs(::Nothing) = (;)
_sigma_kwargs(sigma) = (; sigma)

function evaluate_vrhox!(
        model::KSEModel;
        rho::AbstractArray{<:Real},
        vrho::AbstractArray{<:Real},
        sigma = nothing)
    evaluate_functional!(model.exchange; rho = rho, vrho = vrho, _sigma_kwargs(sigma)...)
end

function evaluate_vrhoc!(
        model::KSEModel;
        rho::AbstractArray{<:Real},
        vrho::AbstractArray{<:Real},
        sigma = nothing)
    evaluate_functional!(model.correlation; rho = rho, vrho = vrho, _sigma_kwargs(sigma)...)
end

function evaluate_vrho!(model::KSEModel;
        rho::AbstractArray{<:Real},
        vrho::AbstractArray{<:Real},
        cache::AbstractArray{<:Real},
        sigma = nothing)
    evaluate_vrhoc!(model; rho = rho, vrho = cache, sigma = sigma)
    evaluate_vrhox!(model; rho = rho, vrho = vrho, sigma = sigma)
    @. vrho += cache
end

function evaluate_zk!(
        model::KSEModel;
        rho::AbstractArray{<:Real},
        zk::AbstractArray{<:Real},
        cache::AbstractArray{<:Real},
        sigma = nothing)
    evaluate_functional!(model.correlation; rho = rho, zk = cache, _sigma_kwargs(sigma)...)
    evaluate_functional!(model.exchange; rho = rho, zk = zk, _sigma_kwargs(sigma)...)
    @. zk += cache
end

"""
    evaluate_vrho_vsigma!(model; rho, sigma, vrho, vsigma, cache_vrho, cache_vsigma)

GGA counterpart of [`evaluate_vrho!`](@ref) that also returns
`vsigma = ∂(ρ·εxc)/∂σ`, needed for the weak-form GGA XC potential matrix (see
`dev/gga/PLAN.md`). Only meaningful when `has_gga(model)`; a `NoFunctional`
component contributes zero to both `vrho` and `vsigma`.
"""
function evaluate_vrho_vsigma!(model::KSEModel;
        rho::AbstractArray{<:Real},
        sigma::AbstractArray{<:Real},
        vrho::AbstractArray{<:Real},
        vsigma::AbstractArray{<:Real},
        cache_vrho::AbstractArray{<:Real},
        cache_vsigma::AbstractArray{<:Real})
    evaluate_functional!(model.correlation; rho = rho, sigma = sigma,
                         vrho = cache_vrho, vsigma = cache_vsigma)
    evaluate_functional!(model.exchange; rho = rho, sigma = sigma,
                         vrho = vrho, vsigma = vsigma)
    @. vrho += cache_vrho
    @. vsigma += cache_vsigma
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

"""
    PBE(; Z::Real, N::Real, nspin::Int = 1)

Convenience constructor for an extended Kohn–Sham model using the PBE GGA
exchange–correlation functional (Libxc `:gga_x_pbe` / `:gga_c_pbe`).

# Arguments
- `Z::Real`: Nuclear charge (e.g., 10.0 for neon).
- `N::Real`: Number of electrons.

# Returns
- `KSEModel` using the PBE GGA exchange and correlation functionals.
"""
function PBE(; Z::Real, N::Real, nspin::Int = 1)
    ex = Functional(:gga_x_pbe, n_spin = nspin)
    ec = Functional(:gga_c_pbe, n_spin = nspin)
    return KSEModel(Z = Z, N = N, ex = ex, ec = ec)
end
