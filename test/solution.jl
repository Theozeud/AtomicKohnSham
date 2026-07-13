# Struct definitions must be at top level (not inside a testset's local scope).
struct _HasFields
    a::Int
    b::String
end

@testset "KSEContext" begin
    mesh = linmesh(0.0, 5.0, 6)
    basis = P1IntLegendreBasis(mesh; ordermax = 3)
    model = KSEModel(Z = 1.0, N = 1.0)
    alg = ODA(tinit = 0.6, aufbau = OptimizedAufbau(max_degen = 1), scftol = 1e-8)
    fim = GaussLegendre(basis, 50)
    ctx = AtomicKohnSham.KSEContext(model, alg, 1, 2, basis, fim)
    @test ctx.model === model
    @test ctx.alg === alg
    @test ctx.lh == 1
    @test ctx.nh == 2
    @test ctx.basis === basis
    @test ctx.fem_integration_method === fim
end

@testset "getfield_or_nothing" begin
    x = _HasFields(3, "hi")
    @test AtomicKohnSham.getfield_or_nothing(x, :a) == 3
    @test AtomicKohnSham.getfield_or_nothing(x, :b) == "hi"
    @test AtomicKohnSham.getfield_or_nothing(x, :nonexistent) === nothing
end

@testset "_occupied_orbitals_summary" begin
    mesh0 = linmesh(0.0, 5.0, 6)
    basis0 = P1IntLegendreBasis(mesh0; ordermax = 3)
    model0 = KSEModel(Z = 1.0, N = 1.0)
    dis0 = KSEDiscretization(basis0, model0; lh = 0, nh = 1)
    @test AtomicKohnSham._occupied_orbitals_summary(dis0, nothing, nothing) === nothing

    mesh = linmesh(0.0, 5.0, 6)
    basis = P1IntLegendreBasis(mesh; ordermax = 3)
    model = KSEModel(Z = 1.0, N = 1.0)
    dis = KSEDiscretization(basis, model; lh = 1, nh = 2)  # 4 orbitals: 1s,2s,2p,3p

    n = zeros(2, 2)     # (lh+1, nh)
    ϵ = zeros(2, 2)
    n[1, 1] = 2.0; ϵ[1, 1] = -0.5    # 1s, occupied
    n[1, 2] = 0.0; ϵ[1, 2] = -0.2    # 2s, empty -> must be excluded
    n[2, 1] = 1.0; ϵ[2, 1] = -0.3    # 2p, occupied
    n[2, 2] = 0.0; ϵ[2, 2] = -0.1    # 3p, empty -> must be excluded

    summary = AtomicKohnSham._occupied_orbitals_summary(dis, n, ϵ)
    @test length(summary) == 2  # only the two nonzero-occupation orbitals
    # sorted by orbital energy ascending: 1s (-0.5) before 2p (-0.3)
    @test summary[1][1] == "1s"
    @test summary[1][2] == -0.5
    @test summary[1][3] == 2.0
    @test summary[2][1] == "2p"
    @test summary[2][2] == -0.3
    @test summary[2][3] == 1.0
end

@testset "eval_nuclear / eval_kinetic_potential (pure formulas)" begin
    # eval_nuclear: -Z/r, with -Inf exactly at r=0
    mesh = linmesh(0.0, 5.0, 6)
    basis = P1IntLegendreBasis(mesh; ordermax = 3)
    model = KSEModel(Z = 3.0, N = 1.0)
    alg = ODA(tinit = 0.6, aufbau = OptimizedAufbau(max_degen = 1), scftol = 1e-8)
    ctx = AtomicKohnSham.KSEContext(model, alg, 0, 1, basis, GaussLegendre(basis, 50))
    dis = KSEDiscretization(basis, model; lh = 0, nh = 1)
    solver = KSESolver(model, dis, alg; maxiter = 1)
    sol = KSESolution(solver, "test")

    X = [0.0, 1.0, 2.0, 5.0]
    Vnuc = eval_nuclear(sol, X)
    @test Vnuc == [-Inf, -3.0, -1.5, -0.6]

    # eval_kinetic_potential(l, X): l(l+1)/(2r^2); 0 at r=0 for l=0, +Inf for l>0
    @test eval_kinetic_potential(0, X) == [0.0, 0.0, 0.0, 0.0]
    Vcent2 = eval_kinetic_potential(2, X)
    @test Vcent2[1] == Inf
    @test Vcent2[2:end] ≈ 2*3 ./ (2 .* X[2:end].^2)
    # the (sol, l, X) method must agree with the pure (l, X) method
    @test eval_kinetic_potential(sol, 1, X) == eval_kinetic_potential(1, X)
end

@testset "eval_orbital/eval_density/eval_hartree/eval_vxc: exact hydrogen cross-check" begin
    # Bare Coulomb (hartree=0, no XC), exactly solvable: the FEM-computed
    # orbital and density must reproduce the analytical 1s solution.
    #   R_10(r) = 2 Z^{3/2} e^{-Zr}      (physical radial wavefunction u(r)/r)
    #   ρ(r)    = R_10(r)^2 / (4π)       (Y_00 = 1/sqrt(4π))
    Z = 1.0
    ex = NoFunctional(n_spin = 1)
    ec = NoFunctional(n_spin = 1)
    model = KSEModel(; Z = Z, N = 1.0, hartree = 0, ex = ex, ec = ec)
    mesh = expmesh(0, 1000, 30; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = 10)
    dis = KSEDiscretization(basis, model; lh = 0, nh = 1)
    aufbau = OptimizedAufbau(max_degen = 1)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-12)
    sol = groundstate(model, dis, alg; maxiter = 200)

    X = [0.5, 1.0, 2.0, 3.0]
    R_exact = 2 .* Z^1.5 .* exp.(-Z .* X)
    R_computed = eval_orbital(sol, 1, 0, X)
    @test isapprox(R_computed, R_exact; rtol = 1e-9)

    rho_exact = R_exact .^ 2 ./ (4π)
    rho_computed = eval_density(sol, X)
    @test isapprox(rho_computed, rho_exact; rtol = 1e-9)

    # eval_orbital via the shell-label string method must agree with the
    # (n, l) method.
    @test eval_orbital(sol, "1s", X) == eval_orbital(sol, 1, 0, X)

    # No exchange-correlation: eval_vxc is identically zero.
    @test all(iszero, eval_vxc(sol, X))

    # hartree=0: W is empty, so eval_hartree must reduce to just the boundary
    # term N/Rmax, not throw or read garbage.
    # Regression test: this used to be an out-of-bounds read of sol.W (masked
    # by @inbounds), only visible with --check-bounds=yes.
    boundary = model.N / last(basis.mesh)
    @test eval_hartree(sol, X) == fill(boundary, length(X))

    # l=0, hartree=0, no XC: the effective potential is exactly the bare
    # nuclear Coulomb potential.
    @test isapprox(eval_effective_potential(sol, 0, X), eval_nuclear(sol, X); rtol = 1e-9)
end

@testset "write_report / Base.show" begin
    ex = NoFunctional(n_spin = 1)
    ec = NoFunctional(n_spin = 1)
    model = KSEModel(; Z = 1, N = 1, ex = ex, ec = ec)
    mesh = expmesh(0, 30, 20; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = 6)
    dis = KSEDiscretization(basis, model; lh = 0, nh = 1)
    alg = ODA(tinit = 0.6, aufbau = OptimizedAufbau(max_degen = 1), scftol = 1e-8)
    sol = groundstate(model, dis, alg; maxiter = 100)

    io = IOBuffer()
    show(io, sol)
    shown = String(take!(io))
    @test occursin("Success", shown)
    @test occursin("Occupation number", shown)
    @test occursin("Sanity checks", shown)
    @test occursin("virial ratio", shown)

    path = tempname()
    @test_nowarn write_report(sol, path)
    report = read(path, String)
    rm(path)
    @test occursin("PARAMETERS", report)
    @test occursin("ENERGIES", report)
    @test occursin("OCCUPIED ORBITALS", report)
    @test occursin("SANITY CHECKS", report)
    @test occursin(string(sol.energies.Etot), report)
end

@testset "Plotting stubs (no CairoMakie loaded)" begin
    # These are catch-all fallbacks: without `using CairoMakie`, they must
    # raise a clear, actionable error rather than a generic MethodError.
    for f in (AtomicKohnSham.plot_density, AtomicKohnSham.plot_orbitals,
        AtomicKohnSham.plot_potentials, AtomicKohnSham.plot_convergence,
        AtomicKohnSham.plot_energy_breakdown)
        err = try
            f()
            nothing
        catch e
            e
        end
        @test err !== nothing
        @test occursin("CairoMakie", sprint(showerror, err))
    end
end
