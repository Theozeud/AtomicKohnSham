"""
    capacity(ch, nspin) -> Int

Number of electrons an orbital of channel `ch` can hold: `2l+1`, times 2 when
the model is spin-unpolarized.
"""
capacity(ch::ChannelData, nspin::Int) = (3 - nspin) * (2ch.l + 1)

"""
    scf_step!(alg::QCFON, solver) -> stopping_criteria

ODA iterations until their stopping criterion falls below `alg.switch_tol`, then
Newton iterations.
"""
function scf_step!(alg::QCFON, solver::KSESolver)
    cache = solver.algcache
    if !cache.newton
        crit = oda_step!(alg.warmup, cache.oda, solver)
        cache.t = cache.oda.t
        solver.niter > 0 && crit < cache.switch_tol && start_newton!(cache, alg, solver)
        return crit
    end
    newton_step!(cache, alg, solver)
end

scf_converged(alg::QCFON, solver::KSESolver) =
    solver.stopping_criteria < alg.scftol || solver.algcache.stalled

function postcomputations!(cache::QCFONCache, solver::KSESolver)
    @unpack discretization, model = solver
    assemble_hartree_pot!(discretization, cache.D; coeff = model.hartree)
end

"""
    start_newton!(cache, alg, solver)

Read the structure of the ODA iterate: the occupied orbitals of each channel
with their occupations, the fractional ones, and the Fermi level (mean energy of
the fractional orbitals, or highest occupied energy if there is none).
"""
function start_newton!(cache::QCFONCache{T}, alg::QCFON, solver::KSESolver) where T
    @unpack discretization = solver
    @unpack nspin = discretization
    @unpack n, ϵ = cache.oda
    empty!(cache.partial)
    ϵocc = T[]
    for (c, ch) in enumerate(cache.channels)
        nk = n[ch.l + 1, :, ch.σ]
        ch.occ = findall(>(0), nk)
        ch.nocc = nk[ch.occ]
        w = capacity(ch, nspin)
        for (j, k) in enumerate(ch.occ)
            push!(ϵocc, ϵ[ch.l + 1, k, ch.σ])
            ch.nocc[j] < w - alg.occ_tol && push!(cache.partial, (c, j))
        end
    end
    cache.μ = isempty(cache.partial) ? maximum(ϵocc) :
        sum(ϵ[cache.channels[c].l + 1, cache.channels[c].occ[j], cache.channels[c].σ]
            for (c, j) in cache.partial) / length(cache.partial)
    for σ in 1:nspin
        @views cache.Dn[:, :, σ] .= _spin(cache.oda.D, σ)
    end
    cache.Dstart .= cache.oda.D
    copy_energies!(cache.Estart, solver.energies)
    cache.newton = true
    cache.t = one(T)
    cache.prevcrit = zero(T)
    @debug "QCFON: Newton starts" occupations = join(["l=$(ch.l) σ=$(ch.σ) n=$(round.(ch.nocc; digits = 6))" for ch in cache.channels], "; ") partial = [(cache.channels[c].l, cache.channels[c].σ,
        cache.channels[c].occ[j], cache.channels[c].nocc[j]) for (c, j) in cache.partial] cache.μ
    nothing
end

"""
    copy_energies!(dest::Energies, src::Energies)

Copy every energy contribution of `src` into `dest`.
"""
function copy_energies!(dest::Energies, src::Energies)
    dest.Etot, dest.Ekin, dest.Ecou = src.Etot, src.Ekin, src.Ecou
    dest.Ehar, dest.Eexc, dest.Ekincor = src.Ehar, src.Eexc, src.Ekincor
    dest
end

"""
    newton_step!(cache, alg, solver) -> stopping_criteria

One Newton iteration: diagonalize the Hamiltonian of the iterate `Dn`, build the
output density matrix with the current occupations, record the corresponding
state and energies, and, unless converged, solve the Newton system by GMRES and
update `Dn`, the fractional occupations and the Fermi level.
"""
function newton_step!(cache::QCFONCache{T}, alg::QCFON, solver::KSESolver) where T
    @unpack discretization, model = solver
    @unpack nspin = discretization
    @unpack Dn, Dout, channels, partial = cache
    assemble_hamiltonian!(discretization, model, nspin == 1 ? _spin(Dn, 1) : Dn)
    diagonalize_channels!(channels, discretization, discretization.ksham.H)
    density_matrix!(Dout, channels)
    while !aufbau_valid(cache, alg, discretization)
        norm(Dout .- Dn) < alg.structure_tol &&
            update_structure!(cache, alg, discretization) && continue
        cache.oda.D .= cache.Dstart
        copy_energies!(solver.energies, cache.Estart)
        cache.newton = false
        cache.switch_tol /= 100
        @debug "QCFON: Aufbau check failed, back to ODA" cache.switch_tol
        return norm(Dout .- Dn)
    end
    partial = cache.partial
    foreach(response_coefficients!, channels)
    density_matrix!(Dout, channels)
    R = Dout .- Dn
    rp = fon_residuals(cache)
    rc = isempty(rp) ? zero(T) : electron_residual(cache, model)
    crit = isempty(rp) ? norm(R) : max(norm(R), norm(rp), abs(rc))
    record_state!(cache, solver)
    crit < alg.scftol && return crit
    if zero(T) < cache.prevcrit < 2crit && crit < sqrt(alg.scftol)
        cache.stalled = true
        @debug "QCFON: residual at its roundoff floor" crit cache.prevcrit
        return crit
    end
    cache.prevcrit = crit

    np = length(rp)
    b = pack_newton(R, rp, rc)
    prepare_newton_operator!(cache, solver)
    tol = max(min(T(0.1), crit) * norm(b), 10eps(T) * norm(b))
    x, info = linsolve(x -> newton_operator(x, cache, solver), b, zero(b),
                       GMRES(krylovdim = alg.krylovdim, maxiter = alg.maxiter_krylov,
                             tol = tol))
    @debug "QCFON Newton step" crit info.converged info.numops info.normres tol
    δD, δλ, δμ = unpack_newton(x, size(Dn), np)
    Dn .+= δD
    for (q, (c, j)) in enumerate(partial)
        np == 0 && break
        channels[c].nocc[j] += capacity(channels[c], nspin) * δλ[q]
    end
    cache.μ += δμ
    crit
end

"""
    fermi_level(cache) -> μ

Fermi level of the current structure: the Newton unknown `cache.μ` when at least
two orbitals are fractional, the energy of the fractional orbital when there is
one, the highest occupied energy otherwise.
"""
function fermi_level(cache::QCFONCache)
    @unpack channels, partial = cache
    length(partial) ≥ 2 && return cache.μ
    length(partial) == 1 && return channels[partial[1][1]].ϵ[channels[partial[1][1]].occ[partial[1][2]]]
    maximum(ch.ϵ[i] for ch in channels for i in ch.occ)
end

"""
    aufbau_valid(cache, alg, discretization) -> Bool

Aufbau check of the current eigenpairs with the current structure, within
`alg.aufbau_tol`: fractional occupations within the capacity of their orbital,
and, among the `nₕ` lowest orbitals of each channel (those the model can occupy,
as for the Aufbau of ODA), empty orbitals not below and full orbitals not above
the Fermi level [`fermi_level`](@ref).
"""
function aufbau_valid(cache::QCFONCache{T}, alg::QCFON,
                      discretization::KSEDiscretization) where T
    @unpack channels, partial = cache
    @unpack nspin, nₕ = discretization
    tol = alg.aufbau_tol
    for (c, j) in partial
        w = capacity(channels[c], nspin)
        if !(-tol ≤ channels[c].nocc[j] ≤ w + tol)
            @debug "QCFON: fractional occupation out of bounds" info = "l=$(channels[c].l) n=$(channels[c].nocc[j])"
            return false
        end
    end
    μ = fermi_level(cache)
    for (c, ch) in enumerate(channels), a in 1:nₕ
        j = findfirst(==(a), ch.occ)
        if j === nothing && ch.ϵ[a] < μ - tol
            @debug "QCFON: empty orbital below the Fermi level" info = "l=$(ch.l) k=$a ϵ=$(ch.ϵ[a]) μ=$μ"
            return false
        elseif j !== nothing && (c, j) ∉ partial && ch.ϵ[a] > μ + tol
            @debug "QCFON: full orbital above the Fermi level" info = "l=$(ch.l) k=$a ϵ=$(ch.ϵ[a]) μ=$μ"
            return false
        end
    end
    true
end

"""
    update_structure!(cache, alg, discretization) -> Bool

Active-set update of the structure (which orbitals are full, fractional, empty),
from the current eigenpairs, when the Aufbau check fails:

1. a fractional occupation below 0 (resp. above the capacity) is set to 0 (resp.
   to the capacity) and the orbital becomes empty (resp. full); the electrons in
   excess are shared among the remaining fractional orbitals, in proportion to
   their capacities, so that the electron number is unchanged;
2. among the `nₕ` lowest orbitals of each channel, an empty orbital below the
   Fermi level becomes fractional (with occupation 0), and a full orbital above it
   becomes fractional (with its full occupation);
3. if this leaves a single fractional orbital that was just added, the highest
   full orbital becomes fractional too, so that electrons can move between them.

The Fermi level is [`fermi_level`](@ref) of the structure before the update, the
same as in [`aufbau_valid`](@ref). Returns `false` when there is no consistent
update (electrons in excess with no fractional orbital left, more than two
orbitals changing status, an update that changes nothing), or after
`alg.max_structure_updates` updates. It is only attempted when the residual is
below `alg.structure_tol` (see [`newton_step!`](@ref)).
"""
function update_structure!(cache::QCFONCache{T}, alg::QCFON,
                           discretization::KSEDiscretization) where T
    @unpack channels = cache
    @unpack nspin, nₕ = discretization
    cache.nupdates >= alg.max_structure_updates && return false
    cache.nupdates += 1
    tol = alg.aufbau_tol
    cap(c) = capacity(channels[c], nspin)
    occ = [Dict(zip(ch.occ, ch.nocc)) for ch in channels]
    frac = Set((c, channels[c].occ[j]) for (c, j) in cache.partial)

    excess = zero(T)
    for (c, k) in collect(frac)
        nk = occ[c][k]
        if nk < 0
            excess += nk
            delete!(occ[c], k)
            delete!(frac, (c, k))
        elseif nk > cap(c)
            excess += nk - cap(c)
            occ[c][k] = cap(c)
            delete!(frac, (c, k))
        end
    end
    if !iszero(excess)
        isempty(frac) && return false
        wsum = sum(cap(c) for (c, _) in frac)
        for (c, k) in frac
            occ[c][k] += excess * cap(c) / wsum
        end
    end

    μ = fermi_level(cache)
    added = Tuple{Int, Int}[]
    for (c, ch) in enumerate(channels), a in 1:nₕ
        if !haskey(occ[c], a) && ch.ϵ[a] < μ - tol
            occ[c][a] = zero(T)
            push!(added, (c, a))
        elseif haskey(occ[c], a) && (c, a) ∉ frac && ch.ϵ[a] > μ + tol
            push!(added, (c, a))
        end
    end
    union!(frac, added)
    if length(frac) == 1 && !isempty(added)
        full = [(c, k) for c in eachindex(channels) for k in keys(occ[c]) if (c, k) ∉ frac]
        isempty(full) && return false
        push!(frac, full[argmax([channels[c].ϵ[k] for (c, k) in full])])
    end

    nchanged = length(added) + length(setdiff(Set((c, channels[c].occ[j]) for (c, j) in cache.partial), frac))
    nchanged > 2 && return false
    unchanged = iszero(excess) && isempty(added) &&
                frac == Set((c, channels[c].occ[j]) for (c, j) in cache.partial)
    unchanged && return false
    empty!(cache.partial)
    for (c, ch) in enumerate(channels)
        ch.occ = sort!(collect(keys(occ[c])))
        ch.nocc = [occ[c][k] for k in ch.occ]
        for (j, k) in enumerate(ch.occ)
            (c, k) ∈ frac && push!(cache.partial, (c, j))
        end
    end
    length(cache.partial) ≥ 2 && length(added) > 0 && (cache.μ = μ)
    cache.prevcrit = zero(T)
    @debug "QCFON: structure updated" fractional = join(["(l=$(channels[c].l), k=$(channels[c].occ[j]), n=$(round(channels[c].nocc[j]; digits = 6)), ϵ=$(round(channels[c].ϵ[channels[c].occ[j]]; digits = 6)))" for (c, j) in cache.partial], " ") excess μ
    true
end

"""
    fon_residuals(cache) -> Vector

Degeneracy residuals `μ - ϵ_p` of the fractional orbitals, when there are at least
two of them (with a single one, its occupation is fixed by the electron number).
"""
function fon_residuals(cache::QCFONCache{T}) where T
    length(cache.partial) < 2 && return T[]
    [cache.μ - cache.channels[c].ϵ[cache.channels[c].occ[j]] for (c, j) in cache.partial]
end

"""
    newton_operator(x, cache, solver) -> y

Matrix-free Newton operator, after [`prepare_newton_operator!`](@ref). With
`x = (δD, δλ, δμ)` (no `δλ, δμ` if fewer than two fractional orbitals):
`y = ((I - χ₀K) δD - Σ_p (2l_p+1) δλ_p u_p u_pᵀ, ⟨u_p|K δD|u_p⟩ - δμ, Σ_p (2l_p+1) δλ_p)`.
"""
function newton_operator(x::AbstractVector, cache::QCFONCache{T}, solver::KSESolver) where T
    @unpack discretization, model = solver
    @unpack Dn, δH, δDχ, channels, partial, ws = cache
    np = length(partial) < 2 ? 0 : length(partial)
    δD, δλ, δμ = unpack_newton(x, size(Dn), np)
    apply_kernel!(δH, ws, discretization, model, Dn, δD; prepare = false)
    apply_chi0!(δDχ, channels, δH)
    yD = δD .- δDχ
    yp = zeros(T, np)
    for (q, (c, j)) in enumerate(partial)
        np == 0 && break
        ch = channels[c]
        @views u = ch.U[:, ch.occ[j]]
        w = capacity(ch, discretization.nspin)
        @views yD[:, :, ch.σ] .-= (w * δλ[q]) .* (u * u')
        @views yp[q] = dot(u, δH[:, :, ch.σ], u) - δμ
    end
    yc = np == 0 ? T[] : [sum(capacity(channels[c], discretization.nspin) * δλ[q]
                              for (q, (c, _)) in enumerate(partial))]
    vcat(vec(yD), yp, yc)
end

"""
    electron_residual(cache, model) -> T

Electron-number residual `N - Σ nocc` of the current occupations. The Newton
steps conserve the electron number only up to the GMRES tolerance, so this
residual is part of the Newton right-hand side, which corrects the drift.
"""
electron_residual(cache::QCFONCache, model::KSEModel) =
    model.N - sum(sum(ch.nocc; init = zero(model.N)) for ch in cache.channels)

"""
    pack_newton(R, rp, rc) -> b

Newton right-hand side: the density-matrix residual `R`, and, with at least two
fractional orbitals, the degeneracy residuals `rp` and the electron-number
residual `rc`.
"""
pack_newton(R::AbstractArray, rp::AbstractVector, rc) =
    isempty(rp) ? vec(copy(R)) : vcat(vec(R), rp, rc)

"""
    unpack_newton(x, dims, np) -> (δD, δλ, δμ)

Split a Newton vector into the density-matrix variation (reshaped to `dims`), the
`np` occupation variations and the Fermi-level variation.
"""
function unpack_newton(x::AbstractVector{T}, dims, np::Int) where T
    m = prod(dims)
    δD = reshape(x[1:m], dims)
    np == 0 && return δD, T[], zero(T)
    δD, x[m+1:m+np], x[m+np+1]
end

"""
    record_state!(cache, solver)

Copy the current orbitals, energies and occupations (lowest `nₕ` of each channel)
and the output density matrix into the shared arrays read by `KSESolution`, and
compute the energies of this state.
"""
function record_state!(cache::QCFONCache{T}, solver::KSESolver) where T
    @unpack discretization, model, energies = solver
    @unpack nₕ, nspin = discretization
    @unpack U, ϵ, n, D, Dout, channels = cache
    fill!(n, zero(T))
    for ch in channels
        @views U[:, :, ch.l + 1, ch.σ] .= ch.U[:, 1:nₕ]
        @views ϵ[ch.l + 1, :, ch.σ] .= ch.ϵ[1:nₕ]
        for (j, k) in enumerate(ch.occ)
            n[ch.l + 1, k, ch.σ] = ch.nocc[j]
        end
    end
    for σ in 1:nspin
        _spin(D, σ) .= @view Dout[:, :, σ]
    end
    energies.Ekin = compute_kinetic_energy(discretization, U, n)
    energies.Ecou = compute_coulomb_energy(discretization, U, n)
    energies.Ehar = compute_hartree_energy(discretization, D)
    energies.Eexc = compute_exc_energy(discretization, model, D)
    energies.Etot = energies.Ekin + energies.Ecou + energies.Ehar + energies.Eexc
    nothing
end

"""
    prepare_newton_operator!(cache, solver)

Precompute what the Newton operator needs and that does not change during one
Newton step: the exchange-correlation kernel of the iterate `Dn`.
"""
function prepare_newton_operator!(cache::QCFONCache, solver::KSESolver)
    @unpack discretization, model = solver
    has_exchcorr(model) && prepare_kernel!(cache.ws, discretization, model, cache.Dn)
    nothing
end
