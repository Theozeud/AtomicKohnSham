using AtomicKohnSham: eval_v2rho2, ResponseWorkspace, apply_kernel!, apply_chi0!,
    empty_channels, diagonalize_channels!, density_matrix!, response_coefficients!,
    assemble_hamiltonian!

function small_lda_problem(Z, N, nspin; lh = 1, nh = 3)
    ex = BuiltinFunctional(:lda_x; nspin = nspin)
    ec = BuiltinFunctional(:lda_c_pw; nspin = nspin)
    model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)
    basis = P1IntLegendreBasis(expmesh(0, 40, 12; s = 1.2); ordermax = 6)
    dis = KSEDiscretization(basis, model; lh = lh, nh = nh,
                            fem_integration_method = GaussLegendre(basis, 200))
    model, dis
end

const ODA_TEST = ODA(tinit = 0.6, scftol = 1e-11,
                     aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5))

@testset "XC kernels against Libxc" begin
    for (lib, rtol2) in ((:lda_x, 1e-12), (:lda_c_pw, 1e-6))
        f1, f2 = BuiltinFunctional(lib; nspin = 1), BuiltinFunctional(lib; nspin = 2)
        x1, x2 = Functional(lib; n_spin = 1), Functional(lib; n_spin = 2)
        for ρ in (1e-4, 1e-2, 1.0, 1e2)
            @test isapprox(eval_v2rho2(f1, ρ),
                           Libxc.evaluate(x1; rho = [ρ], derivatives = 2).v2rho2[1];
                           rtol = 1e-12)
            for ζ in (-0.5, 0.0, 0.8)
                ρu, ρd = ρ * (1 + ζ) / 2, ρ * (1 - ζ) / 2
                lx = Libxc.evaluate(x2; rho = reshape([ρu, ρd], 2, 1),
                                    derivatives = 2).v2rho2[:, 1]
                k = collect(eval_v2rho2(f2, ρu, ρd))
                @test norm(k - lx) ≤ rtol2 * norm(lx)
            end
        end
    end
end

@testset "Linear response against finite differences" begin
    model, dis = small_lda_problem(10, 10, 1)
    sol = groundstate(model, dis, ODA_TEST; maxiter = 300)
    model2, dis2 = small_lda_problem(10, 9.7, 1)
    sol2 = groundstate(model2, dis2, ODA_TEST; maxiter = 300)
    nb = dis.Nₕ
    D = reshape(copy(sol.D), nb, nb, 1)
    δD = reshape(sol2.D .- sol.D, nb, nb, 1)
    δH = zeros(nb, nb, 1)
    apply_kernel!(δH, ResponseWorkspace(dis), dis, model, D, δD)
    h = 1e-4
    meanfield(X) = (assemble_hamiltonian!(dis, model, X);
                    dis.ksham.H[:, :, 1, 1] .- dis.ksham.Hfix[1])
    fd = (meanfield(sol.D .+ h .* δD[:, :, 1]) .- meanfield(sol.D .- h .* δD[:, :, 1])) ./ 2h
    @test norm(δH[:, :, 1] - fd) ≤ 1e-5 * norm(fd)

    assemble_hamiltonian!(dis, model, sol.D)
    H = copy(dis.ksham.H)
    channels = diagonalize_channels!(empty_channels(dis), dis, H)
    for ch in channels
        nk = sol.n[ch.l + 1, :, ch.σ]
        ch.occ = findall(>(0), nk)
        ch.nocc = nk[ch.occ]
        response_coefficients!(ch)
    end
    δDχ = zeros(nb, nb, 1)
    apply_chi0!(δDχ, channels, δH)
    hh = 1e-5
    Dp, Dm = zeros(nb, nb, 1), zeros(nb, nb, 1)
    density_matrix!(Dp, diagonalize_channels!(deepcopy(channels), dis, H, hh .* δH))
    density_matrix!(Dm, diagonalize_channels!(deepcopy(channels), dis, H, -hh .* δH))
    @test norm(δDχ - (Dp .- Dm) ./ 2hh) ≤ 1e-4 * norm(δDχ)
end

@testset "QCFON reproduces ODA" begin
    for (Z, N, nspin) in ((10, 10, 1), (6, 6, 2), (17, 17.5, 2))
        model, dis = small_lda_problem(Z, N, nspin)
        so = groundstate(model, dis, ODA_TEST; maxiter = 300)
        model, dis = small_lda_problem(Z, N, nspin)
        sq = groundstate(model, dis, QCFON(scftol = 1e-11); maxiter = 100)
        @test sq.success == "SUCCESS"
        @test sq.niter < so.niter
        @test abs(sq.energies.Etot - so.energies.Etot) < 1e-9
        @test sq.occupied == so.occupied || all(isapprox(a[2], b[2]; atol = 1e-8) &&
                                                isapprox(a[3], b[3]; atol = 1e-8)
                                                for (a, b) in zip(sq.occupied, so.occupied))
    end
end
