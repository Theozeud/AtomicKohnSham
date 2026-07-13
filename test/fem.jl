@testset "Mesh" begin
    # Every mesh generator must honor its (a, b, n, T) contract: right element
    # type, endpoints pinned exactly, and the requested point count.
    a, b, n, T, s = 0.0, 1.0, 10, Float64, 0.9

    # --- Linear mesh ---
    @test_nowarn linmesh(a, b, n; T = T)
    lin_mesh = linmesh(a, b, n; T = T)
    @test eltype(lin_mesh) == T
    @test first(lin_mesh) == a
    @test last(lin_mesh) == b
    @test length(lin_mesh) == n

    # --- Geometric mesh ---
    @test_nowarn geometricmesh(a, b, n; T = T, s = s)
    geo_mesh = geometricmesh(a, b, n; T = T, s = s)
    @test eltype(geo_mesh) == T
    @test first(geo_mesh) == a
    @test last(geo_mesh) == b
    @test length(geo_mesh) == n

    # --- Polynomial mesh ---
    @test_nowarn polynomialmesh(a, b, n; T = T, s = s)
    poly_mesh = polynomialmesh(a, b, n; T = T, s = s)
    @test eltype(poly_mesh) == T
    @test first(poly_mesh) == a
    @test last(poly_mesh) == b
    @test length(poly_mesh) == n

    # --- Exponential mesh ---
    @test_nowarn expmesh(a, b, n; T = T, s = s)
    exp_mesh = expmesh(a, b, n; T = T, s = s)
    @test eltype(exp_mesh) == T
    @test first(exp_mesh) == a
    @test last(exp_mesh) == b
    @test length(exp_mesh) == n

    # --- API ---
    @test_nowarn print(linmesh)
end

@testset "Mesh Robustness" begin
    # --- findindex: cell localisation at nodes, midpoints, and out-of-domain points ---
    lin = linmesh(0.0, 5.0, 6)  # points 0,1,2,3,4,5 -> 5 cells
    @test AtomicKohnSham.findindex(lin, 0.0) == 1  # left boundary -> first cell
    @test AtomicKohnSham.findindex(lin, 5.0) == 5  # right boundary -> last cell
    @test AtomicKohnSham.findindex(lin, 2.0) == 3  # interior node -> its right-owning cell
    for k in 1:5
        mid = (lin[k] + lin[k + 1]) / 2
        @test AtomicKohnSham.findindex(lin, mid) == k
    end
    # Out-of-domain points must return an out-of-range sentinel rather than throw
    # (callers are responsible for guarding on it; see the "FEM Basis" out-of-domain tests).
    @test_nowarn AtomicKohnSham.findindex(lin, -1.0)
    @test AtomicKohnSham.findindex(lin, -1.0) == 0
    @test AtomicKohnSham.findindex(lin, 6.0) == length(lin) + 1
    # Floating-point roundoff just below the lower bound (the scenario that used
    # to crash evaluate!) must resolve to the same "below range" sentinel, not a
    # valid-looking-but-wrong cell index.
    @test AtomicKohnSham.findindex(lin, -1e-12) == 0

    # --- Mesh generation must produce strictly increasing points, with both
    #     endpoints pinned exactly, across generator types, grading exponents,
    #     and resolutions (no zero-width or reversed cells). ---
    a, b = 0.0, 5.0
    configs = [
        (linmesh, NamedTuple()),
        (geometricmesh, (s = 0.85,)),
        (geometricmesh, (s = 1.15,)),
        (polynomialmesh, (s = 1.5,)),
        (polynomialmesh, (s = 0.7,)),
        (expmesh, (s = 1.7,)),
        (expmesh, (s = 1.2,)),
    ]
    for (meshgen, kwargs) in configs, n in (5, 20, 100)
        m = meshgen(a, b, n; T = Float64, kwargs...)
        @test issorted(m.points)
        @test all(diff(m.points) .> 0)
        @test first(m) == a
        @test last(m) == b
        @test length(m) == n
    end
end

@testset "FEM Basis" begin
    mesh = linmesh(0.0, 2.0, 5)  # points 0, 0.5, 1.0, 1.5, 2.0 -> 4 cells
    basis = P1IntLegendreBasis(mesh; ordermax = 3)
    n = length(basis)

    # size accounting: (ordermax-1) bubbles per cell + one shared hat per interior node
    @test n == (3 - 1) * 4 + (4 - 1)

    xs = [0.25, 0.5, 0.75, 1.2, 1.999]

    # --- Scalar evaluation pb(i, x) vs vector evaluation pb(I, x): two independent
    #     code paths (the scalar method never touches the evaluate! cache) that must agree. ---
    for x in xs
        @test basis(1:n, x) == [basis(i, x) for i in 1:n]
    end

    # --- Whole-basis evaluation pb(x) must match pb(1:length(pb), x).
    #     Regression test: pb(x) previously threw UndefVarError (referenced an undefined `I`). ---
    for x in xs
        @test basis(x) == basis(1:n, x)
    end

    # --- Reconstructing a single basis function from its one-hot coefficient vector
    #     through evaluate!/evaluate (the coeffs+X vectorized path) must match the
    #     ground-truth single-function evaluation pb(i, x).
    #     `evaluate` is qualified below because Libxc also exports a function named
    #     `evaluate`; test/runtests.jl loads both packages, so the bare name is
    #     ambiguous there even though it resolves fine in isolation. ---
    X = [0.1, 0.6, 1.1, 1.6]
    for i in 1:n
        coeffs = zeros(n)
        coeffs[i] = 1.0
        @test AtomicKohnSham.evaluate(basis, coeffs, X) ≈ [basis(i, x) for x in X]
    end

    # --- Continuity across a cell boundary (interior mesh node x=1.0, shared by two
    #     cells): left and right limits of a nontrivial coefficient combination must agree. ---
    coeffs = ones(n)
    left = AtomicKohnSham.evaluate(basis, coeffs, [1.0 - 1e-6])[1]
    right = AtomicKohnSham.evaluate(basis, coeffs, [1.0 + 1e-6])[1]
    at_node = AtomicKohnSham.evaluate(basis, coeffs, [1.0])[1]
    @test isapprox(left, at_node; atol = 1e-3)
    @test isapprox(right, at_node; atol = 1e-3)

    # --- Out-of-domain points: basis functions have compact support, so points
    #     outside [first(mesh), last(mesh)] must evaluate to zero, not throw.
    #     Regression test: evaluate/evaluate! previously crashed with a KeyError/BoundsError
    #     for such points (e.g. tiny floating-point noise just below r=0). ---
    below = first(mesh) - 1.0
    above = last(mesh) + 1.0
    @test basis(1, below) == 0
    @test basis(1, above) == 0
    @test_nowarn AtomicKohnSham.evaluate(basis, ones(n), [below, above])
    @test AtomicKohnSham.evaluate(basis, ones(n), [below, above]) == [0.0, 0.0]

    # The specific failure mode that motivated the fix: a point an ulp below the
    # mesh's lower bound (plausible floating-point roundoff at r=0).
    just_below = first(mesh) - 1e-12
    @test_nowarn AtomicKohnSham.evaluate(basis, ones(n), [just_below])
    @test AtomicKohnSham.evaluate(basis, ones(n), [just_below]) == [0.0]
end

@testset "FEM Matrices" begin
    # Basis under test: must match test/generate_fem_reference_data.jl exactly.
    # Rerun that script (not part of the test suite) to refresh the reference
    # file below after a deliberate change to basis/integration code.
    a, b, n, T, s = 0.0, 1.0, 3, Float64, 0.9
    mesh = polynomialmesh(a, b, n; T = T, s = s)
    basis = P1IntLegendreBasis(mesh, T; ordermax = 2)

    ref = deserialize(joinpath(@__DIR__, "reference_data", "fem_matrices.jls"))
    tol = 1e-12

    # --- Overlap (mass) matrix, plain and weighted by 1/x, 1/x^2 ---
    @test isapprox(mass_matrix(basis), ref["M0"]; atol = tol)
    @test isapprox(Matrix(sparse_mass_matrix(basis)), ref["M0"]; atol = tol)
    @test isapprox(mass_matrix(basis, -1), ref["M1"]; atol = tol)
    @test isapprox(Matrix(sparse_mass_matrix(basis, -1)), ref["M1"]; atol = tol)
    @test isapprox(mass_matrix(basis, -2), ref["M2"]; atol = tol)
    @test isapprox(Matrix(sparse_mass_matrix(basis, -2)), ref["M2"]; atol = tol)

    # --- Stiffness matrix ---
    @test isapprox(stiffness_matrix(basis), ref["A"]; atol = tol)
    @test isapprox(Matrix(sparse_stiffness_matrix(basis)), ref["A"]; atol = tol)

    # --- Mass tensor, plain and weighted ---
    @test isapprox(mass_tensor(basis), ref["T0"]; atol = tol)
    @test isapprox(mass_tensor(basis, -1), ref["T1"]; atol = tol)
    @test isapprox(mass_tensor(basis, -2), ref["T2"]; atol = tol)
end

@testset "FEM Matrices Properties" begin
    # Structural/analytical checks that hold for any mesh/order, independent of
    # the fixed reference values in "FEM Matrices" above (which only pin down
    # one specific setup).
    mesh = expmesh(0.0, 3.0, 6; T = Float64, s = 1.3)  # 5 cells
    basis = P1IntLegendreBasis(mesh; ordermax = 4)
    n = length(basis)

    # --- Mass matrix: symmetric, and positive definite (Gram matrix of
    #     independent basis functions), plain and weighted by 1/x, 1/x^2 ---
    M = mass_matrix(basis)
    @test isapprox(M, M'; atol = 1e-10)
    @test minimum(eigvals(Symmetric(M))) > 0

    M1 = mass_matrix(basis, -1)
    @test isapprox(M1, M1'; atol = 1e-8)
    @test minimum(eigvals(Symmetric(M1))) > 0

    M2 = mass_matrix(basis, -2)
    @test isapprox(M2, M2'; atol = 1e-6)
    @test minimum(eigvals(Symmetric(M2))) > 0

    # --- Stiffness matrix: symmetric, and strictly positive definite (no
    #     boundary hats -> no representable constants -> no null space) ---
    A = stiffness_matrix(basis)
    @test isapprox(A, A'; atol = 1e-10)
    @test minimum(eigvals(Symmetric(A))) > 0

    # --- Mass tensor T[i,j,k] = ∫ φᵢφⱼφₖ must be symmetric under any
    #     permutation of its three indices (the integrand is a product,
    #     order doesn't matter). ---
    Tt = mass_tensor(basis)
    for (i, j, k) in [(1, 2, 3), (4, 7, 2), (n, 1, n - 1), (3, 3, 5), (n, n, n)]
        v = Tt[i, j, k]
        @test Tt[j, i, k] ≈ v
        @test Tt[i, k, j] ≈ v
        @test Tt[k, j, i] ≈ v
        @test Tt[j, k, i] ≈ v
        @test Tt[k, i, j] ≈ v
    end

    # --- Cross-check the GaussLegendre/FunWeight quadrature assembly path
    #     against ExactIntegration: a constant weight of 1 must reproduce the
    #     plain mass matrix. This exercises the same code path used for
    #     Hartree/XC assembly in the SCF loop, validated against the
    #     independently-computed exact result. ---
    const_weight = AtomicKohnSham.FunWeight(
        X -> ones(eltype(X), length(X)); is_vectorized = true, is_inplace = false)
    Agl = zeros(eltype(basis), n, n)
    AtomicKohnSham.fill_mass_matrix!(
        basis, Agl; weight = const_weight, method = GaussLegendre(basis, 60))
    @test isapprox(Agl, M; atol = 1e-8, rtol = 1e-8)
end
