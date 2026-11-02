# Stages 2-3: linear response, on lithium (LSDA, Slater + PW92) and neon (LDA).
#   - apply_kernel!: Hamiltonian variation K δD, against central finite
#     differences of assemble_hamiltonian! (Hartree + exchange-correlation part);
#   - apply_chi0!: density-matrix variation χ₀ δH at fixed occupations, against
#     central finite differences of the density matrix of H ± h δH.
# Double64 throughout, so that finite differences are not limited by roundoff
# (in Float64 they stall at 1e-6 through the diagonalization).
# Run: julia stage2_response.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", "examples"))
using AtomicKohnSham, LinearAlgebra, Printf, DoubleFloats
const T = Double64
using AtomicKohnSham: ResponseWorkspace, apply_kernel!, apply_chi0!, empty_channels,
    diagonalize_channels!, density_matrix!, response_coefficients!, assemble_hamiltonian!

"""Converged ODA ground state on a small discretization."""
function setup(Z, N, nspin; lh = 1, nh = 3)
    ex = BuiltinFunctional(:lda_x; nspin = nspin)
    ec = BuiltinFunctional(:lda_c_pw; nspin = nspin)
    model = KSEModel(; Z = T(Z), N = T(N), ex = ex, ec = ec)
    basis = P1IntLegendreBasis(expmesh(zero(T), T(40), 12; s = T(1.2)), T; ordermax = 6)
    dis = KSEDiscretization(basis, model; lh = lh, nh = nh, eigensolver = FullEigenSolver(),
                            fem_integration_method = GaussLegendre(basis, 200))
    alg = ODA(tinit = T(0.6), aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5), scftol = T(1e-20))
    sol = groundstate(model, dis, alg; maxiter = 300)
    model, dis, sol
end

as3(D, nspin) = reshape(copy(D), size(D, 1), size(D, 2), nspin)

"""Hartree + exchange-correlation part of the Hamiltonian at D, per spin."""
function mean_field(dis, model, D, nspin)
    assemble_hamiltonian!(dis, model, D)
    H = dis.ksham.H
    cat([H[:, :, 1, σ] .- dis.ksham.Hfix[1] for σ in 1:nspin]...; dims = 3)
end

relerr(a, b) = Float64(norm(a - b) / norm(b))

for (name, Z, N, nspin) in (("Li LSDA", 3, 3, 2), ("Ne LDA", 10, 10, 1))
    model, dis, sol = setup(Z, N, nspin)
    nb = dis.Nₕ
    D = as3(sol.D, nspin)
    _, _, sol2 = setup(Z, N - 1, nspin)
    δD = as3(sol2.D, nspin) .- D
    ws = ResponseWorkspace(dis)

    δH = zeros(T, nb, nb, nspin)
    apply_kernel!(δH, ws, dis, model, D, δD)
    h = T(1e-8)
    Dp, Dm = D .+ h .* δD, D .- h .* δD
    arg(X) = nspin == 1 ? X[:, :, 1] : X
    fd = (mean_field(dis, model, arg(Dp), nspin) .- mean_field(dis, model, arg(Dm), nspin)) ./ 2h
    @printf("%-8s stage 2  K δD vs FD:   rel err %.1e\n", name, relerr(δH, fd))

    assemble_hamiltonian!(dis, model, arg(D))
    H = copy(dis.ksham.H)
    channels = empty_channels(dis)
    diagonalize_channels!(channels, dis, H)
    for ch in channels
        nk = sol.n[ch.l + 1, :, ch.σ]
        ch.occ = findall(>(0), nk)
        ch.nocc = nk[ch.occ]
        response_coefficients!(ch)
    end
    δDχ = zeros(T, nb, nb, nspin)
    apply_chi0!(δDχ, channels, δH)
    hh = T(1e-8)
    Dpp, Dmm = zeros(T, nb, nb, nspin), zeros(T, nb, nb, nspin)
    density_matrix!(Dpp, diagonalize_channels!(deepcopy(channels), dis, H, hh .* δH))
    density_matrix!(Dmm, diagonalize_channels!(deepcopy(channels), dis, H, -hh .* δH))
    @printf("%-8s stage 3  χ₀ δH vs FD:  rel err %.1e\n", name, relerr(δDχ, (Dpp .- Dmm) ./ 2hh))
end
