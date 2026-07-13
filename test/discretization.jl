@testset "Discretization Indexing" begin
    # Bare discretizations (no init_cache! / SCF needed): convert_index,
    # convert_index_nl, and degeneracy only read lₕ/nₕ/nspin off the struct.
    mesh = linmesh(0.0, 5.0, 6)
    basis = P1IntLegendreBasis(mesh; ordermax = 3)

    for (nspin, lh, nh) in [(1, 2, 3), (2, 2, 3), (1, 0, 1), (2, 3, 2)]
        model = KSEModel(
            Z = 1.0, N = 1.0, ex = NoFunctional(nspin), ec = NoFunctional(nspin))
        dis = KSEDiscretization(basis, model; lh = lh, nh = nh)
        total = (lh + 1) * nh * nspin

        # --- convert_index is a bijection onto its valid (l,k,σ) range ---
        seen = Set{Tuple{Int, Int, Int}}()
        for idx in 1:total
            l, k, σ = AtomicKohnSham.convert_index(dis, idx)
            @test 0 <= l <= lh
            @test 1 <= k <= nh
            @test 1 <= σ <= nspin
            @test (l, k, σ) ∉ seen  # never repeats
            push!(seen, (l, k, σ))
        end
        @test length(seen) == total  # covers every valid combination exactly once

        # --- convert_index_nl agrees with convert_index (n = k+l, same l/σ) ---
        for idx in 1:total
            l, k, σ = AtomicKohnSham.convert_index(dis, idx)
            nl = AtomicKohnSham.convert_index_nl(dis, idx)
            if nspin == 1
                @test nl == (k + l, l)
            else
                @test nl == (k + l, l, σ)
            end
        end

        # --- degeneracy formula: 4l+2 (unpolarized) or 2l+1 (polarized) ---
        for idx in 1:total
            l, _, _ = AtomicKohnSham.convert_index(dis, idx)
            expected = nspin == 1 ? 4l + 2 : 2l + 1
            @test AtomicKohnSham.degeneracy(dis, idx) == expected
        end

        # --- independent closed-form check: total capacity across all l (and,
        #     for nspin=2, both spin channels) for a single radial index k is
        #     2(lh+1)^2 -- the standard "full shell" electron-capacity formula
        #     (sum of the first (lh+1) odd numbers is (lh+1)^2, times 2 for spin
        #     however it's tracked), derived analytically rather than by
        #     re-summing degeneracy() itself. ---
        cap_per_k = sum(idx -> AtomicKohnSham.convert_index(dis, idx)[2] == 1 ?
                               AtomicKohnSham.degeneracy(dis, idx) : 0, 1:total)
        @test cap_per_k == 2 * (lh + 1)^2
    end
end
