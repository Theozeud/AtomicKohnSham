@testset "Mesh" begin
    a = 0.0
    b = 1.0
    n = 10
    T = Float64
    s = 0.9

    # LINEAR MESH
    @test_nowarn linmesh(a, b, n; T = T)
    lin_mesh = linmesh(a, b, n; T = T)
    @test eltype(lin_mesh) == T
    @test first(lin_mesh) == a
    @test last(lin_mesh) == b
    @test length(lin_mesh) == n

    # GEOMETRIC MESH
    @test_nowarn geometricmesh(a, b, n; T = T, s = s)
    geo_mesh = geometricmesh(a, b, n; T = T, s = s)
    @test eltype(geo_mesh) == T
    @test first(geo_mesh) == a
    @test last(geo_mesh) == b
    @test length(geo_mesh) == n

    # POLYNOMIAL MESH
    @test_nowarn polynomialmesh(a, b, n; T = T, s = s)
    poly_mesh = polynomialmesh(a, b, n; T = T, s = s)
    @test eltype(poly_mesh) == T
    @test first(poly_mesh) == a
    @test last(poly_mesh) == b
    @test length(poly_mesh) == n

    # EXPONENTIAL MESH
    @test_nowarn polynomialmesh(a, b, n; T = T, s = s)
    exp_mesh = polynomialmesh(a, b, n; T = T, s = s)
    @test eltype(exp_mesh) == T
    @test first(exp_mesh) == a
    @test last(exp_mesh) == b
    @test length(exp_mesh) == n

    # API
    @test_nowarn print(linmesh)
end

@testset "FEM Basis" begin end

@testset "FEM Matrices" begin
    a = 0.0
    b = 1.0
    n = 3
    T = Float64
    s = 0.9
    exp_mesh = polynomialmesh(a, b, n; T = T, s = s)

    basis = P1IntLegendreGenerator(exp_mesh, T; ordermax = 2)

    # Overlap Matrix
    @test_nowarn mass_matrix(basis)
    @test_nowarn sparse_mass_matrix(basis)
    M0 = mass_matrix(basis)
    M0_sparse = sparse_mass_matrix(basis)

    @test norm(M0 - [1.333333333333333 -0.17862891042271553 -0.1547044229106178;
                -0.17862891042271553 0.07145156416908623 0.0;
                -0.1547044229106178 0.0 0.06188176916424714]) < 1e-14

    @test norm(M0_sparse - sparse([1, 2, 3, 1, 2, 1, 3],
        [1, 1, 1, 2, 2, 3, 3],
        [1.333333333333333, -0.17862891042271553, -0.1547044229106178,
            -0.17862891042271553, 0.07145156416908623, -0.1547044229106178,
            0.06188176916424714], 3, 3)) < 1e-14

    # Weight Mass Matrix
    @test_nowarn mass_matrix(basis, -1)
    @test_nowarn sparse_mass_matrix(basis, -1)
    M1 = mass_matrix(basis, -1)
    M1_sparse = sparse_mass_matrix(basis, -1)

    @test norm(M1 - [2.965986233780117 -0.6666666666666666 -0.21796077205123415;
                -0.6666666666666666 0.3333333333333333 0.0;
                -0.21796077205123415 0.0 0.081665704871926]) < 1e-14

    @test norm(M1_sparse - sparse([1, 2, 3, 1, 2, 1, 3],
        [1, 1, 1, 2, 2, 3, 3],
        [2.965986233780117, -0.6666666666666666, -0.21796077205123415,
            -0.6666666666666666, 0.3333333333333333, -0.21796077205123415,
            0.081665704871926], 3, 3)) < 1e-14

    @test_nowarn mass_matrix(basis, -2)
    @test_nowarn sparse_mass_matrix(basis, -2)
    M2 = mass_matrix(basis, -2)
    M2_sparse = sparse_mass_matrix(basis, -2)

    @test norm(M2 - [8.996555397028684 -3.7321319661472296 -0.312103917626415;
                -3.7321319661472296 2.488087977431486 0.0;
                -0.312103917626415 0.0 0.10925872461476452]) < 1e-14

    @test norm(M2_sparse - sparse([1, 2, 3, 1, 2, 1, 3],
        [1, 1, 1, 2, 2, 3, 3],
        [8.996555397028684, -3.7321319661472296, -0.312103917626415,
            -3.7321319661472296, 2.488087977431486, -0.312103917626415,
            0.10925872461476452], 3, 3)) < 1e-14

    # Stiffness Matrix
    @test_nowarn stiffness_matrix(basis)
    @test_nowarn sparse_stiffness_matrix(basis)
    A = stiffness_matrix(basis)
    A_sparse = sparse_stiffness_matrix(basis)

    @test norm(A - [16.082849673076296 0.0 0.0;
                0.0 2.488087977431486 0.0;
                0.0 0.0 2.8728619135939453]) < 1e-14

    @test norm(A_sparse - sparse([1, 2, 3],
        [1, 2, 3],
        [16.082849673076296, 2.488087977431486, 2.8728619135939453],
        3, 3)) < 1e-14

    # Mass Tensor
    @test_nowarn mass_tensor(basis)
    T0 = mass_tensor(basis)

    @test norm(T0 - [2.0 -0.21435469250725864 -0.18564530749274138;
                -0.21435469250725864 0.07145156416908623 0.0;
                -0.18564530749274138 0.0 0.06188176916424714;;;
                -0.21435469250725864 0.07145156416908623 0.0;
                0.07145156416908623 -0.03062209892960838 0.0;
                0.0 0.0 0.0;;;
                -0.18564530749274138 0.0 0.06188176916424714;
                0.0 0.0 0.0;
                0.06188176916424714 0.0 -0.02652075821324877]) < 1e-14

    @test_nowarn mass_tensor(basis, -1)
    T1 = mass_tensor(basis, -1)

    @test norm(T1 - [4.162717590124432 -0.6666666666666666 -0.27259013435861623;
                -0.6666666666666666 0.2666666666666667 0.0;
                -0.27259013435861623 0.0 0.08525477309337263;;;
                -0.6666666666666666 0.2666666666666667 0.0;
                0.2666666666666667 -0.13333333333333336 0.0;
                0.0 0.0 0.0;;;
                -0.27259013435861623 0.0 0.08525477309337263;
                0.0 0.0 0.0;
                0.08525477309337263 0.0 -0.03489421349767301]) < 1e-14

    @test_nowarn mass_tensor(basis, -2)
    T2 = mass_tensor(basis, -2)

    @test norm(T2 - [9.904639026510079 -2.4880879774314866 -0.40569038602330093;
                -2.4880879774314866 1.2440439887157433 0.0;
                -0.40569038602330093 0.0 0.1189064032503698;;;
                -2.4880879774314866 1.2440439887157433 0.0;
                1.2440439887157433 -0.7464263932294459 0.0;
                0.0 0.0 0.0;;;
                -0.40569038602330093 0.0 0.1189064032503698;
                0.0 0.0 0.0;
                0.1189064032503698 0.0 -0.046399038294079875]) < 1e-14
end
