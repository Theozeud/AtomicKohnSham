@testset "Slater Exchange" begin
    func1 = AtomicKohnSham.SlaterXα(1)  # spin-unpolarized
    func2 = AtomicKohnSham.SlaterXα(2)  # spin-polarized

    # --- exact closed-form values (Dirac/Slater exchange, Cx = 3/4*(3/pi)^(1/3)) ---
    Cx = 3/4 * (3/π)^(1/3)
    ρ = 0.42
    @test isapprox(AtomicKohnSham.eval_zk(func1, ρ), -Cx * ρ^(1/3); atol = 1e-14)
    @test isapprox(AtomicKohnSham.eval_vrho(func1, ρ), -(3/π)^(1/3) * ρ^(1/3); atol = 1e-14)

    # --- eval_vrho is the derivative of ρ*eval_zk(ρ) w.r.t. ρ (finite-difference check) ---
    h = 1e-6
    d_num = (ρ * AtomicKohnSham.eval_zk(func1, ρ + h) -
             (ρ - h) * AtomicKohnSham.eval_zk(func1, ρ - h)) / (2h)
    # centered around ρ using ρ+h and ρ-h with matching prefactors
    energy_density(x) = x * AtomicKohnSham.eval_zk(func1, x)
    d_num2 = (energy_density(ρ + h) - energy_density(ρ - h)) / (2h)
    @test isapprox(d_num2, AtomicKohnSham.eval_vrho(func1, ρ); atol = 1e-6)

    # --- spin-polarized formulas reduce to the unpolarized ones at ρup=ρdown=ρ/2
    #     (the defining self-consistency property of the spin-scaling relation) ---
    for ρtest in [0.1, 0.42, 3.7]
        @test isapprox(AtomicKohnSham.eval_zk(func2, ρtest/2, ρtest/2),
                        AtomicKohnSham.eval_zk(func1, ρtest); atol = 1e-13)
        @test isapprox(AtomicKohnSham.eval_vrho_up(func2, ρtest/2, ρtest/2),
                        AtomicKohnSham.eval_vrho(func1, ρtest); atol = 1e-13)
        @test isapprox(AtomicKohnSham.eval_vrho_down(func2, ρtest/2, ρtest/2),
                        AtomicKohnSham.eval_vrho(func1, ρtest); atol = 1e-13)
    end

    # --- fully spin-polarized limit (ρdown=0): reduces to the unpolarized
    #     functional evaluated at 2ρup (a fully spin-polarized system behaves
    #     like a spin-unpolarized system at twice the density, per spin channel) ---
    ρup = 0.8
    @test isapprox(AtomicKohnSham.eval_vrho_up(func2, ρup, 0.0),
                    AtomicKohnSham.eval_vrho(func1, 2ρup); atol = 1e-13)
    @test AtomicKohnSham.eval_vrho_down(func2, ρup, 0.0) == AtomicKohnSham.eval_vrho(func1, 0.0)

    # --- vacuum edge case: no NaN/Inf from the ρup+ρdown=0 division in eval_zk ---
    @test AtomicKohnSham.eval_zk(func2, 0.0, 0.0) == 0.0
    @test isfinite(AtomicKohnSham.eval_vrho_up(func2, 0.0, 0.0))

    # --- known symmetry: eval_zk is symmetric under swapping ρup <-> ρdown ---
    @test AtomicKohnSham.eval_zk(func2, 0.3, 0.7) == AtomicKohnSham.eval_zk(func2, 0.7, 0.3)
end

@testset "BuiltinFunctional Dispatch" begin
    # --- constructor validates nspin and no longer throws UndefVarError ---
    @test_nowarn BuiltinFunctional(:lda_x; nspin = 1)
    @test_nowarn BuiltinFunctional(:lda_x; nspin = 2)
    @test_throws Exception BuiltinFunctional(:lda_x; nspin = 3)
    @test_throws Exception BuiltinFunctional(:unknown_functional; nspin = 1)

    f1 = BuiltinFunctional(:lda_x; nspin = 1)
    @test f1 isa AtomicKohnSham.SlaterXα
    @test f1.n_spin == 1
    @test AtomicKohnSham.is_lda(f1)

    # --- NoFunctional always evaluates to zero, regardless of density ---
    nf = NoFunctional(n_spin = 1)
    rho = [0.1, 0.5, 2.0]
    zk = zeros(3)
    vrho = zeros(3)
    AtomicKohnSham.evaluate_functional!(nf; rho = rho, zk = zk, vrho = vrho)
    @test all(iszero, zk)
    @test all(iszero, vrho)

    # --- evaluate_functional!/evaluate_functional on SlaterXα (nspin=1) match
    #     the raw eval_zk/eval_vrho formulas pointwise ---
    func = BuiltinFunctional(:lda_x; nspin = 1)
    zk2 = zeros(3)
    vrho2 = zeros(3)
    AtomicKohnSham.evaluate_functional!(func; rho = rho, zk = zk2, vrho = vrho2)
    @test zk2 ≈ [AtomicKohnSham.eval_zk(func, r) for r in rho]
    @test vrho2 ≈ [AtomicKohnSham.eval_vrho(func, r) for r in rho]

    # --- spin-polarized (nspin=2) evaluation no longer throws (previously
    #     eval_zk/eval_vrho_up/eval_vrho_down were undefined for this path) ---
    funcpol = BuiltinFunctional(:lda_x; nspin = 2)
    rho2 = [0.2 0.3 0.1; 0.1 0.2 0.4]  # 2 spin channels x 3 points
    zk3 = zeros(3)
    vrho3 = zeros(2, 3)
    @test_nowarn AtomicKohnSham.evaluate_functional!(funcpol; rho = rho2, zk = zk3, vrho = vrho3)
    @test zk3 ≈ [AtomicKohnSham.eval_zk(funcpol, rho2[1, i], rho2[2, i]) for i in 1:3]
    @test vrho3[1, :] ≈ [AtomicKohnSham.eval_vrho_up(funcpol, rho2[1, i], rho2[2, i]) for i in 1:3]
    @test vrho3[2, :] ≈ [AtomicKohnSham.eval_vrho_down(funcpol, rho2[1, i], rho2[2, i]) for i in 1:3]
end

@testset "KSEModel" begin
    model_none = KSEModel(Z = 1.0, N = 1.0)
    @test !AtomicKohnSham.has_exchange(model_none)
    @test !AtomicKohnSham.has_correlation(model_none)
    @test !AtomicKohnSham.has_exchcorr(model_none)
    @test model_none.nspin == 1

    model_ex = Slater(Z = 1.0, N = 1.0)
    @test AtomicKohnSham.has_exchange(model_ex)
    @test !AtomicKohnSham.has_correlation(model_ex)
    @test AtomicKohnSham.has_exchcorr(model_ex)

    model_rhf = RHF(Z = 2.0, N = 2.0)
    @test !AtomicKohnSham.has_exchcorr(model_rhf)
    @test model_rhf.hartree == 0 || model_rhf.hartree == 1  # documented convention

    # --- nspin propagates from the exchange/correlation functionals to the model ---
    model_pol = Slater(Z = 1.0, N = 1.0; nspin = 2)
    @test model_pol.nspin == 2
end
