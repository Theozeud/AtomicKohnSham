"""
    QCFON(; scftol, warmup = ODA(...), switch_tol = 1e-2, krylovdim = 30,
            maxiter_krylov = 10, occ_tol = 1e-10, aufbau_tol = 1e-6,
            max_structure_updates = 20, structure_tol = 1e-3)

Quadratically convergent SCF for fractional occupation numbers (Cancès, Kudin,
Scuseria, Turinici, J. Chem. Phys. 118, 2003), in its radial form: a Newton method
on the density matrix, `(I - χ₀K) δD = D_out[D] - D`, where `χ₀` is the
independent-particle response computed by a sum over states in each channel
`(l, σ)` and `K` the Hartree-exchange-correlation kernel. When several channels
carry a fractional orbital, their occupations and the Fermi level are additional
unknowns (degenerate levels stay degenerate, the electron number is fixed).

The structure of the solution (fully and partially occupied orbitals) is taken from
`warmup` iterations, run until their stopping criterion falls below `switch_tol`.
Each Newton system is solved by GMRES (`KrylovKit.linsolve`), without assembling
any matrix. LDA and LSDA only.

# Keywords
- `scftol`: tolerance on the Newton residual `‖D_out[D] - D‖` (and on the
  degeneracy residuals of fractional levels). Since this is a true residual, not
  a step length, it has a roundoff floor that grows with the atom; the iterations
  also stop when the residual is below `√scftol` and no longer decreases (by less
  than a factor 2 in one step), i.e. when it has reached that floor.
- `warmup`: first-order algorithm (ODA) used before the Newton iterations.
- `switch_tol`: stopping criterion of `warmup` below which Newton starts.
- `krylovdim`, `maxiter_krylov`: GMRES Krylov dimension and number of restarts.
- `occ_tol`: tolerance to classify an occupation as full.
- `aufbau_tol`: tolerance of the Aufbau check of the Newton iterates. When it
  fails (an empty orbital below the Fermi level, a full one above it, or a
  fractional occupation leaving its range), the structure is updated as in an
  active-set method (see [`update_structure!`](@ref)); if no consistent update
  exists, or after `max_structure_updates` updates, the iterations go back to
  `warmup`, with a switch tolerance divided by 100.
- `max_structure_updates`: maximal number of structure updates (guards against
  cycling).
- `structure_tol`: the structure is updated only when the current residual is
  below `structure_tol` and at most two orbitals change status (Newton is a local
  method: far from the solution, the iterations go back to `warmup` instead).
"""
struct QCFON{T <: Real, OA <: ODA} <: SCFAlgorithm
    scftol::T
    warmup::OA
    switch_tol::T
    krylovdim::Int
    maxiter_krylov::Int
    occ_tol::T
    aufbau_tol::T
    max_structure_updates::Int
    structure_tol::T
    function QCFON(; scftol::Real,
                   warmup::ODA = ODA(tinit = 0.6, scftol = scftol,
                                     aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5)),
                   switch_tol::Real = 1e-2, krylovdim::Int = 30, maxiter_krylov::Int = 10,
                   occ_tol::Real = 1e-10, aufbau_tol::Real = 1e-6,
                   max_structure_updates::Int = 20, structure_tol::Real = 1e-3)
        T = typeof(scftol)
        new{T, typeof(warmup)}(scftol, warmup, T(switch_tol), krylovdim, maxiter_krylov,
                               T(occ_tol), T(aufbau_tol), max_structure_updates,
                               T(structure_tol))
    end
end

"""
    QCFONCache

Workspace of [`QCFON`](@ref): the ODA workspace used for the first iterations
(whose `D`, `U`, `ϵ`, `n` arrays are shared, so that a solution reads them as for
ODA), the Newton iterate `Dn`, the full eigendecomposition of each channel, the
fractional orbitals `partial` (channel index, position among the occupied
orbitals of that channel), the Fermi level `μ`, the response buffers, the
current switch tolerance (reduced after a failed Aufbau check), and the previous
Newton residual with the roundoff-floor flag, the number of structure updates,
and the ODA iterate and its energies at the start of Newton (restored when going
back to ODA, whose line search relies on the energies of the previous iterate).
"""
mutable struct QCFONCache{T <: Real, OC <: ODACache, DT, UT, ET, NT, WS} <: SCFCache
    oda::OC
    D::DT
    U::UT
    ϵ::ET
    n::NT
    newton::Bool
    Dn::Array{T, 3}
    Dout::Array{T, 3}
    δH::Array{T, 3}
    δDχ::Array{T, 3}
    channels::Vector{ChannelData{T}}
    partial::Vector{Tuple{Int, Int}}
    μ::T
    ws::WS
    t::T
    switch_tol::T
    prevcrit::T
    stalled::Bool
    nupdates::Int
    Dstart::DT
    Estart::Energies{T}
end

function create_cache_alg(alg::QCFON, discretization::KSEDiscretization{T},
                          model::KSEModel) where T
    oda = create_cache_alg(alg.warmup, discretization, model)
    @unpack Nₕ, nspin = discretization
    z() = zeros(T, Nₕ, Nₕ, nspin)
    ws = ResponseWorkspace(discretization)
    QCFONCache{T, typeof(oda), typeof(oda.D), typeof(oda.U), typeof(oda.ϵ), typeof(oda.n),
               typeof(ws)}(oda, oda.D, oda.U, oda.ϵ, oda.n, false, z(), z(), z(), z(),
                           empty_channels(discretization), Tuple{Int, Int}[], zero(T), ws,
                           zero(T), T(alg.switch_tol), zero(T), false, 0, copy(oda.D),
                           Energies(T))
end
