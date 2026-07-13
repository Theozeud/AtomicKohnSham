@testset "Line Search Energy" begin
    # --- degenerate case: both endpoints have identical energies (a=b=0), so
    #     E(t) is constant and any t in [0,1] is a minimizer. ---
    tmin, Emin = AtomicKohnSham.line_search_energy(
        1.0, 1.0, 2.0, 2.0, 0.5, 0.5, 0.5; nl = false)
    @test 0 <= tmin <= 1
    @test isapprox(Emin, 1.0 + 2.0 + 0.5; atol = 1e-12)  # A + H at either endpoint

    # --- pure quadratic with a known interior minimum: a=2, b=-3, c=2 -> t*=3/4 ---
    energy_kin0, energy_kin1 = 0.0, 0.0
    energy_cou0, energy_cou1 = 0.0, 1.0
    energy_har0, energy_har1 = 1.0, 1.0
    energy_har01 = 0.0
    tmin, Emin = AtomicKohnSham.line_search_energy(
        energy_kin0, energy_kin1, energy_cou0, energy_cou1,
        energy_har0, energy_har1, energy_har01; nl = false)
    @test isapprox(tmin, 3/4; atol = 1e-12)
    @test isapprox(Emin, 2*tmin^2 - 3*tmin + 2; atol = 1e-12)

    # --- a<=0 (concave/linear): minimum sits at whichever endpoint is lower,
    #     never in the interior. ---
    tmin_edge, _ = AtomicKohnSham.line_search_energy(
        0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0; nl = false)
    @test tmin_edge == 0.0 || tmin_edge == 1.0

    # --- nl (Brent) path with F ≡ 0 must reproduce the quadratic-only result ---
    tmin_nl, Emin_nl = AtomicKohnSham.line_search_energy(
        energy_kin0, energy_kin1, energy_cou0, energy_cou1,
        energy_har0, energy_har1, energy_har01, t -> 0.0; nl = true)
    @test isapprox(tmin_nl, tmin; atol = 1e-6)
    @test isapprox(Emin_nl, Emin; atol = 1e-6)
end

# Shared bare discretizations (no init_cache!/SCF needed for occupation-assignment tests).
function _bare_discretization(; nspin = 1, lh = 2, nh = 3, N = 4.0)
    mesh = linmesh(0.0, 5.0, 6)
    basis = P1IntLegendreBasis(mesh; ordermax = 3)
    model = KSEModel(Z = 1.0, N = N, ex = NoFunctional(nspin), ec = NoFunctional(nspin))
    dis = KSEDiscretization(basis, model; lh = lh, nh = nh)
    model, dis
end

@testset "FrozenAufbau" begin
    model, dis = _bare_discretization(; nspin = 1, lh = 1, nh = 2, N = 4.0)

    af = FrozenAufbau(Dict("1s" => 2.0, "2s" => 1.5))
    cache = AtomicKohnSham.create_cache_aufbau(dis, model, af)
    n = AtomicKohnSham.zero_occupation_numbers(dis)
    ϵ = AtomicKohnSham.zero_orbitals_energies(dis)  # ϵ is ignored by FrozenAufbau
    U = zeros(1, 1, 1, 1)
    AtomicKohnSham.aufbau!(n, ϵ, U, model, dis, cache, 0)

    l1, k1 = 0, 1  # "1s" -> n=1,l=0 -> k=n-l=1
    l2, k2 = 0, 2  # "2s" -> n=2,l=0 -> k=n-l=2
    @test n[l1 + 1, k1] == 2.0
    @test n[l2 + 1, k2] == 1.5
    @test cache.postcomputations == false

    # --- occupations are copied verbatim regardless of orbital energies ---
    ϵ2 = AtomicKohnSham.zero_orbitals_energies(dis)
    ϵ2 .= 100.0  # wildly different energies
    n2 = AtomicKohnSham.zero_occupation_numbers(dis)
    AtomicKohnSham.aufbau!(n2, ϵ2, U, model, dis, cache, 5)
    @test n2 == n

    # --- spin-suffix mismatch is now a hard error (was previously a silent
    #     @error followed by a confusing BoundsError deeper in the code) ---
    model_pol, dis_pol = _bare_discretization(; nspin = 2, lh = 0, nh = 1, N = 1.0)
    af_nospin = FrozenAufbau(Dict("1s" => 1.0))
    @test_throws Exception AtomicKohnSham.create_cache_aufbau(dis_pol, model_pol, af_nospin)

    af_spin_on_unpolarized = FrozenAufbau(Dict("1sUP" => 1.0))
    @test_throws Exception AtomicKohnSham.create_cache_aufbau(dis, model, af_spin_on_unpolarized)
end

@testset "OptimizedAufbau" begin
    # --- pure helper: collect_degenerate_block! groups consecutive equal-energy
    #     orbitals up to max_degen, and stops at the first energy gap ---
    ϵ = [0.0, 0.0, 1.0, 1.0, 1.0]
    indices_sort = [1, 2, 3, 4, 5]  # already sorted by ϵ
    indices_block = zeros(Int, 2)
    nb = AtomicKohnSham.collect_degenerate_block!(indices_block, indices_sort, ϵ, 1, 1e-6, 2)
    @test nb == 2
    @test indices_block == [1, 2]
    nb2 = AtomicKohnSham.collect_degenerate_block!(indices_block, indices_sort, ϵ, 3, 1e-6, 2)
    @test nb2 == 2  # capped at max_degen even though 3 orbitals share ϵ=1.0
    @test indices_block == [3, 4]

    # --- fill_full_block!/fill_partial_single! write exactly what's asked ---
    n = zeros(5)
    AtomicKohnSham.fill_full_block!(n, [2, 6], [1, 2], 2)
    @test n[1] == 2 && n[2] == 6
    AtomicKohnSham.fill_partial_single!(n, 3, 1.25)
    @test n[3] == 1.25

    # --- full aufbau! on a spin-unpolarized system: exact electron-count
    #     conservation, occupations within [0, degeneracy], lowest orbitals
    #     filled first ---
    model, dis = _bare_discretization(; nspin = 1, lh = 1, nh = 2, N = 5.0)
    aufbau = OptimizedAufbau()
    cache = AtomicKohnSham.create_cache_aufbau(dis, model, aufbau)
    n = AtomicKohnSham.zero_occupation_numbers(dis)
    ϵ = AtomicKohnSham.zero_orbitals_energies(dis)
    ϵ .= reshape([-2.0, -1.0, 0.0, 1.0], size(ϵ))  # (l+1,k) = (1,1),(2,1),(1,2),(2,2)
    U = zeros(1, 1, 1, 1)
    AtomicKohnSham.aufbau!(n, ϵ, U, model, dis, cache, 0)
    @test isapprox(sum(n), 5.0; atol = 1e-10)
    @test all(0 .<= n .<= 6)  # max degeneracy here is 4*1+2=6 for l=1
    # lowest-energy orbital (l=0,k=1, index (1,1), ϵ=-2.0) must be fully occupied
    @test n[1, 1] == 2  # degeneracy(l=0)=2, and it's the lowest energy: fully filled

    # --- degeneracy > max_degen is explicitly rejected, not silently mishandled ---
    ϵ3 = AtomicKohnSham.zero_orbitals_energies(dis)
    fill!(ϵ3, 0.0)  # every orbital exactly degenerate -> triggers the 3-fold case
    n3 = AtomicKohnSham.zero_occupation_numbers(dis)
    model3, dis3 = _bare_discretization(; nspin = 1, lh = 2, nh = 1, N = 1.0)
    aufbau3 = OptimizedAufbau(max_degen = 3)
    cache3 = AtomicKohnSham.create_cache_aufbau(dis3, model3, aufbau3)
    ϵ3b = AtomicKohnSham.zero_orbitals_energies(dis3)  # all zero -> all 3 orbitals degenerate
    n3b = AtomicKohnSham.zero_occupation_numbers(dis3)
    U3 = zeros(1, 1, 1, 1)
    @test_throws Exception AtomicKohnSham.aufbau!(n3b, ϵ3b, U3, model3, dis3, cache3, 0)
end

@testset "SmearedAufbau" begin
    # --- constructor validates Temp > 0 ---
    @test_throws Exception SmearedAufbau(Temp = 0.0)
    @test_throws Exception SmearedAufbau(Temp = -0.1)
    @test_nowarn SmearedAufbau(Temp = 0.01)

    model, dis = _bare_discretization(; nspin = 1, lh = 1, nh = 2, N = 3.0)
    aufbau = SmearedAufbau(Temp = 0.05)
    cache = AtomicKohnSham.create_cache_aufbau(dis, model, aufbau)
    n = AtomicKohnSham.zero_occupation_numbers(dis)
    ϵ = AtomicKohnSham.zero_orbitals_energies(dis)
    ϵ .= reshape([-2.0, -1.0, 0.0, 1.0], size(ϵ))
    U = zeros(1, 1, 1, 1)
    AtomicKohnSham.aufbau!(n, ϵ, U, model, dis, cache, 0)

    # --- electron count is conserved exactly (up to the bisection tolerance) ---
    @test isapprox(sum(n), 3.0; atol = 1e-8)
    # --- every occupation lies within [0, degeneracy] ---
    for i in eachindex(ϵ)
        l, _, _ = AtomicKohnSham.convert_index(dis, i)
        @test 0 <= n[i] <= AtomicKohnSham.degeneracy(dis, i)
    end
    # --- occupation is monotonically non-increasing in orbital energy
    #     (Fermi-Dirac is a monotonic function of ε at fixed μ, Temp) ---
    ϵflat = vec(ϵ)
    nflat = vec(n)
    perm = sortperm(ϵflat)
    fractional_occ = nflat[perm] ./ [AtomicKohnSham.degeneracy(dis, i) for i in perm]
    @test issorted(fractional_occ; rev = true)

    # --- larger Temp spreads occupations further from a sharp step (higher
    #     entropy / less step-like fractional occupation on the frontier orbital) ---
    aufbau_hot = SmearedAufbau(Temp = 1.0)
    cache_hot = AtomicKohnSham.create_cache_aufbau(dis, model, aufbau_hot)
    n_hot = AtomicKohnSham.zero_occupation_numbers(dis)
    AtomicKohnSham.aufbau!(n_hot, ϵ, U, model, dis, cache_hot, 0)
    @test isapprox(sum(n_hot), 3.0; atol = 1e-8)
    # the highest-energy orbital should pick up more (nonzero) occupation at
    # higher Temp than at low Temp, since smearing is stronger
    highest_idx = last(sortperm(ϵflat))
    @test n_hot[highest_idx] > n[highest_idx]
end
