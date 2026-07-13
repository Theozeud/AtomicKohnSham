@testset "LogBook" begin
    lb = LogBook(Float64)
    @test lb.stopping_criteria == Float64[]
    @test lb.Etot == Float64[]
    @test eltype(lb.stopping_criteria) == Float64
end

# A minimal mock callback: counts how many times it was invoked, without
# needing a real KSESolver (callback! only needs *a* value with the fields it
# unpacks, and this mock doesn't unpack anything). Defined at top level since
# Julia doesn't allow struct definitions inside a local (testset) scope.
mutable struct _CountingCallback
    count::Int
end
AtomicKohnSham.callback!(cb::_CountingCallback, ::AtomicKohnSham.KSESolver) = (cb.count += 1; nothing)

@testset "CallbackSet" begin
    @test CallbackSet().callbacks == ()

    cb1 = _CountingCallback(0)
    cb2 = _CountingCallback(0)
    cbset = CallbackSet((cb1, cb2))
    # Exercise through a real solver so callback! sees a genuine KSESolver.
    mesh = linmesh(0.0, 5.0, 6)
    basis = P1IntLegendreBasis(mesh; ordermax = 3)
    model = KSEModel(Z = 1.0, N = 1.0)
    dis = KSEDiscretization(basis, model; lh = 0, nh = 1)
    alg = ODA(tinit = 0.6, aufbau = OptimizedAufbau(max_degen = 1), scftol = 1e-8)

    # callback!(::Nothing, ...) is a no-op, given a genuine solver
    solver = KSESolver(model, dis, alg; maxiter = 10)
    @test_nowarn AtomicKohnSham.callback!(nothing, solver)

    @test_nowarn groundstate(model, dis, alg; maxiter = 10, callback = cbset)
    @test cb1.count > 0
    @test cb1.count == cb2.count  # both callbacks fire once per iteration, in lockstep
end

@testset "LogFileCallback / write_log_header" begin
    mesh = expmesh(0, 30, 20; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = 6)
    model = KSEModel(Z = 1.0, N = 1.0)
    dis = KSEDiscretization(basis, model; lh = 0, nh = 1)
    alg = ODA(tinit = 0.6, aufbau = OptimizedAufbau(max_degen = 1), scftol = 1e-8)

    io = IOBuffer()
    AtomicKohnSham.init_cache!(dis, model)
    write_log_header(io, model, alg, dis)
    header = String(take!(io))
    @test occursin("PARAMETERS", header)
    @test occursin("ITERATIONS", header)
    @test occursin("Z (nuclear charge)", header)

    io2 = IOBuffer()
    sol = groundstate(model, dis, alg; maxiter = 50,
        callback = CallbackSet((LogFileCallback(io2),)))
    log_lines = split(String(take!(io2)), '\n'; keepempty = false)
    @test length(log_lines) == sol.niter  # one line per SCF iteration
    @test all(l -> occursin("iter =", l) && occursin("Etot =", l), log_lines)
    # the recorded Etot on the last line must match the final energy
    @test occursin(string(sol.energies.Etot), log_lines[end])
end

@testset "KSESolver / groundstate status and naming" begin
    mesh = expmesh(0, 30, 20; s = 1.2)
    basis = P1IntLegendreBasis(mesh; ordermax = 6)
    aufbau = OptimizedAufbau(max_degen = 1)

    # --- fresh KSESolver starts at niter=0 with zeroed state ---
    model = KSEModel(Z = 1.0, N = 1.0)
    dis = KSEDiscretization(basis, model; lh = 0, nh = 1)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-8)
    solver = KSESolver(model, dis, alg; maxiter = 100)
    @test solver.niter == 0
    @test solver.stopping_criteria == 0.0
    @test solver.maxiter == 100
    @test solver.energies.Etot == 0.0

    # --- SUCCESS vs MAXITERS status, driven purely by whether the SCF loop
    #     reaches scftol before running out of iterations ---
    dis_conv = KSEDiscretization(basis, model; lh = 0, nh = 1)
    sol_ok = groundstate(model, dis_conv, alg; maxiter = 100)
    @test sol_ok.success == "SUCCESS"
    @test sol_ok.niter < 100

    dis_stuck = KSEDiscretization(basis, model; lh = 0, nh = 1)
    alg_impossible = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-15)
    sol_stuck = groundstate(model, dis_stuck, alg_impossible; maxiter = 1)
    @test sol_stuck.success == "MAXITERS"
    @test sol_stuck.niter == 1

    # --- naming: neutral atom, cation, anion (integer Z), and the
    #     non-integer-Z fallback (regression test: this used to print the
    #     literal, non-interpolated string "Z=Z, N=N") ---
    for (Z, N, expect) in [
        (2.0, 2.0, "Helium"),
        (2.0, 1.0, "Helium1.0+"),
        (1.0, 2.0, "Hydrogen1.0-"),
    ]
        m = KSEModel(Z = Z, N = N)
        d = KSEDiscretization(basis, m; lh = 0, nh = 1)
        sol = groundstate(m, d, alg; maxiter = 100)
        @test sol.name == expect
    end

    model_frac = KSEModel(Z = 1.5, N = 1.0)
    dis_frac = KSEDiscretization(basis, model_frac; lh = 0, nh = 1)
    sol_frac = groundstate(model_frac, dis_frac, alg; maxiter = 100)
    @test sol_frac.name == "Z=1.5, N=1.0"
    @test !occursin("Z=Z", sol_frac.name)  # the literal-string bug, if it ever came back
end

@testset "Energies" begin
    e = Energies(Float64)
    @test e.Etot == 0.0
    @test e.Ekin == 0.0
    @test e.Ecou == 0.0
    @test e.Ehar == 0.0
    @test e.Eexc == 0.0
    @test e.Ekincor == 0.0
    @test eltype(e.Etot) == Float64
end
