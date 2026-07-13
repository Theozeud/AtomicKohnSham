@testset "Orbital Labels" begin
    # --- parse_shell / shell_string are inverses ---
    for s in ["1s", "2s", "2p", "3d", "10d", "7i"]
        n, l = AtomicKohnSham.parse_shell(s)
        @test AtomicKohnSham.shell_string((n, l)) == s
    end
    for (s, σ) in [("1sUP", 1), ("1sDOWN", 2), ("1s↑", 1), ("1s↓", 2)]
        n, l, σparsed = AtomicKohnSham.parse_shell(s)
        @test σparsed == σ
    end

    # --- round-trip through shell_string in both styles ---
    for t in [(2, 0), (3, 1), (4, 2)]
        @test AtomicKohnSham.parse_shell(AtomicKohnSham.shell_string(t)) == t
    end
    for t in [(2, 0, 1), (2, 0, 2), (4, 2, 1)]
        @test AtomicKohnSham.parse_shell(AtomicKohnSham.shell_string(t; style = :word)) == t
        word_str = AtomicKohnSham.shell_string(t; style = :word)
        @test occursin("UP", word_str) || occursin("DOWN", word_str)
    end

    # --- l labels s,p,d,f,g,h,i map to 0..6 ---
    for (label, l) in zip("spdfghi", 0:6)
        @test AtomicKohnSham.parse_shell("$(l + 1)$label") == (l + 1, l)
    end

    # --- edge cases: n<=l must error, as must unknown labels/spin suffixes ---
    @test_throws Exception AtomicKohnSham.parse_shell("1p")  # n=1 <= l=1
    @test_throws Exception AtomicKohnSham.parse_shell("2z")  # unknown label
    @test_throws Exception AtomicKohnSham.parse_shell("1sLEFT")  # unknown spin suffix
    @test_throws Exception AtomicKohnSham.shell_string((1, 1))  # n<=l
    @test_throws Exception AtomicKohnSham.shell_string((1, 2, 3))  # bad tuple length
end

@testset "Periodic Tables" begin
    # --- both lookup tables cover exactly Z=1..118, with no gaps ---
    @test Set(keys(AtomicKohnSham.ATOMIC_NUMBER_TO_NAME)) == Set(1:118)
    @test Set(keys(AtomicKohnSham.ATOMIC_NUMBER_TO_SYMBOL)) == Set(1:118)

    # --- names and symbols are unique (no duplicate entries) ---
    @test length(Set(values(AtomicKohnSham.ATOMIC_NUMBER_TO_NAME))) == 118
    @test length(Set(values(AtomicKohnSham.ATOMIC_NUMBER_TO_SYMBOL))) == 118

    # --- spot-check well-known elements ---
    @test AtomicKohnSham.ATOMIC_NUMBER_TO_NAME[1] == "Hydrogen"
    @test AtomicKohnSham.ATOMIC_NUMBER_TO_SYMBOL[1] == "H"
    @test AtomicKohnSham.ATOMIC_NUMBER_TO_NAME[2] == "Helium"
    @test AtomicKohnSham.ATOMIC_NUMBER_TO_SYMBOL[2] == "He"
    @test AtomicKohnSham.ATOMIC_NUMBER_TO_NAME[21] == "Scandium"
    @test AtomicKohnSham.ATOMIC_NUMBER_TO_SYMBOL[21] == "Sc"
    @test AtomicKohnSham.ATOMIC_NUMBER_TO_NAME[118] == "Oganesson"
end

@testset "Free Allocation" begin
    # --- flexible_zeros collapses the "extra" dimension exactly when it's 1 ---
    @test size(AtomicKohnSham.flexible_zeros(Float64, (3, 4), 1)) == (3, 4)
    @test size(AtomicKohnSham.flexible_zeros(Float64, (3, 4), 2)) == (3, 4, 2)
    @test size(AtomicKohnSham.flexible_zeros(Float64, 1, (3, 4))) == (3, 4)
    @test size(AtomicKohnSham.flexible_zeros(Float64, 2, (3, 4))) == (2, 3, 4)
    @test all(iszero, AtomicKohnSham.flexible_zeros(Float64, (3, 4), 2))

    # --- tensor_matrix_dict!/tensor_vector_dict!: hand-built small tensor,
    #     checked against direct summation over the dict entries ---
    F = Dict((1, 1, 1) => 2.0, (1, 2, 1) => 3.0, (2, 1, 2) => 5.0, (2, 2, 2) => 7.0)
    D = [1.0 2.0; 3.0 4.0]

    B = zeros(2)
    AtomicKohnSham.tensor_matrix_dict!(B, D, F)
    expected_B = zeros(2)
    for ((i, j, m), val) in F
        expected_B[m] += D[i, j] * val
    end
    @test B == expected_B

    # DUP/DDOWN variant sums the two densities before contracting
    DUP = D
    DDOWN = [0.5 0.0; 0.0 0.5]
    Bspin = zeros(2)
    AtomicKohnSham.tensor_matrix_dict!(Bspin, DUP, DDOWN, F)
    @test isapprox(Bspin, B .+ [F[(1,1,1)]*DDOWN[1,1] + F[(1,2,1)]*DDOWN[2,1],
                                 F[(2,1,2)]*DDOWN[1,2] + F[(2,2,2)]*DDOWN[2,2]])

    Dvec = [10.0, 20.0]
    Bmat = zeros(2, 2)
    AtomicKohnSham.tensor_vector_dict!(Bmat, Dvec, F)
    expected_Bmat = zeros(2, 2)
    for ((i, j, m), val) in F
        expected_Bmat[i, j] += Dvec[m] * val
    end
    @test Bmat == expected_Bmat

    # tensor_matrix_dict! and tensor_vector_dict! are transposes of the same
    # contraction: dot(B_from_matrix_dict, Dvec) == dot(D, B_from_vector_dict)
    @test dot(B, Dvec) ≈ dot(D, Bmat)
end

@testset "Sparse Tools" begin
    # --- symmetrize_sparse! averages (i,j)/(j,i) pairs, leaves the diagonal alone ---
    A = sparse([1, 2, 1, 3], [2, 1, 1, 3], [1.0, 3.0, 5.0, 9.0], 3, 3)
    # A[1,2]=1.0, A[2,1]=3.0 (asymmetric off-diagonal pair), A[1,1]=5.0, A[3,3]=9.0
    AtomicKohnSham.symmetrize_sparse!(A)
    @test A[1, 2] == A[2, 1] == 2.0  # average of 1.0 and 3.0
    @test A[1, 1] == 5.0  # diagonal untouched
    @test A[3, 3] == 9.0  # untouched entry with no symmetric partner
    @test issymmetric(Matrix(A))

    @test_throws Exception AtomicKohnSham.symmetrize_sparse!(sparse([1, 2], [1, 1], [1.0, 2.0], 2, 1))

    # --- scale_sparse! multiplies only the stored (nonzero) entries ---
    B = sparse([1, 2], [1, 2], [2.0, 4.0], 2, 2)
    AtomicKohnSham.scale_sparse!(B, 3.0)
    @test Matrix(B) == [6.0 0.0; 0.0 12.0]
end
