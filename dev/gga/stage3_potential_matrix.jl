# Stage 3 (see PLAN.md): validates the assembled GGA XC potential matrix
# (assemble_exc!'s VxcUP/VxcDOWN, via the weak-form derivation
#     Vxc = M[vrho - 2K/r] + Mmix[K] + Mmix[K]ᵀ, K = 2vσ·ρ'
# ) against a direct analytic-vs-numerical-derivative check: for a fixed,
# arbitrary (not necessarily self-consistent) density matrix D,
#     VxcUP[i,j] == (compute_exc_energy(D + h·e_ij) - compute_exc_energy(D - h·e_ij)) / 2h
# This is the required gate before trusting the GGA assembly (PLAN.md Stage 3):
# it catches any sign/factor mistake in the ∂Exc/∂D_ij formula directly,
# without laundering it through an SCF loop.
#
# Run: julia --project=. dev/gga/stage3_potential_matrix.jl

using AtomicKohnSham
using Test

# ---------------------------------------------------------------------
# nspin = 1 (PBE on Hydrogen): reuse an LDA-converged density matrix as an
# arbitrary, genuinely non-trivial D (same trick as stage2_density_gradient.jl)
# ---------------------------------------------------------------------
Z, N = 1.0, 1.0
Rmax, Nmesh, ordermax = 12.0, 8, 6
mesh = expmesh(0, Rmax, Nmesh; s = 1.2)
basis = P1IntLegendreBasis(mesh; ordermax = ordermax)

ex_lda = BuiltinFunctional(:lda_x; nspin = 1)
ec_lda = NoFunctional(1)
model_lda = KSEModel(; Z = Z, N = N, ex = ex_lda, ec = ec_lda)
dis_lda = KSEDiscretization(basis, model_lda; lh = 0, nh = 1)
aufbau = OptimizedAufbau(max_degen = 1)
alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = 1e-10)
sol = groundstate(model_lda, dis_lda, alg; maxiter = 100)
D = sol.D

model = PBE(Z = Z, N = N, nspin = 1)
dis = KSEDiscretization(basis, model; lh = 0, nh = 1)

AtomicKohnSham.assemble_exc!(dis, model, D)
VxcUP = collect(dis.ksham.VxcUP)

n = length(basis)
h = 1e-6
pairs = [(1, 1), (2, 2), (3, 3), (1, 2), (2, 3), (1, 4), (4, 6), (6, 6), (2, 8), (8, 8)]
pairs = filter(((i, j),) -> i <= n && j <= n, pairs)

@testset "GGA VxcUP[i,j] vs finite difference of compute_exc_energy (nspin=1, PBE)" begin
    for (i, j) in pairs
        # eval_density_gradient!/optimized_eval_density_gradient! (Stage 2)
        # deliberately exploit D's symmetry (q' = 2Σ D_ij φᵢ'φⱼ, valid only
        # for symmetric D -- see PLAN.md's derivation and stage2's docstring).
        # A single-entry, asymmetry-breaking perturbation (D[i,j] alone) would
        # silently feed that shortcut an unphysical D and give a wrong ρ' --
        # so the perturbation must keep D symmetric: bump both D[i,j] and
        # D[j,i] together. For i != j that changes the bilinear form through
        # both D_ij and D_ji, so the derivative doubles: compare against
        # 2*VxcUP[i,j] (VxcUP symmetric) instead of VxcUP[i,j] directly.
        Dp = copy(D)
        Dm = copy(D)
        Dp[i, j] += h; i != j && (Dp[j, i] += h)
        Dm[i, j] -= h; i != j && (Dm[j, i] -= h)
        Ep = AtomicKohnSham.compute_exc_energy(dis, model, Dp)
        Em = AtomicKohnSham.compute_exc_energy(dis, model, Dm)
        fd = (Ep - Em) / (2h)
        expected = i == j ? VxcUP[i, i] : 2 * VxcUP[i, j]
        @test isapprox(expected, fd; atol = 1e-6, rtol = 1e-4)
    end
end

println("Sample check: VxcUP[1,2]=$(VxcUP[1,2])")

# ---------------------------------------------------------------------
# nspin = 2 (PBE on Lithium, spin-polarized): same gate, but now also
# exercising the cross spin term K↑ = 2vσ↑↑·ρ↑' + vσ↑↓·ρ↓' (and ↑↔↓).
# ---------------------------------------------------------------------
Z2, N2 = 3.0, 3.0
ex_lsda = BuiltinFunctional(:lda_x; nspin = 2)
ec_lsda = BuiltinFunctional(:lda_c_pw; nspin = 2)
model_lsda = KSEModel(; Z = Z2, N = N2, ex = ex_lsda, ec = ec_lsda)
dis_lsda = KSEDiscretization(basis, model_lsda; lh = 0, nh = 2)
aufbau2 = OptimizedAufbau(max_degen = 1)
alg2 = ODA(tinit = 0.6, aufbau = aufbau2, scftol = 1e-10)
sol2 = groundstate(model_lsda, dis_lsda, alg2; maxiter = 100)
D2 = sol2.D

model2 = PBE(Z = Z2, N = N2, nspin = 2)
dis2 = KSEDiscretization(basis, model2; lh = 0, nh = 2)
AtomicKohnSham.assemble_exc!(dis2, model2, D2)
VxcUP2 = collect(dis2.ksham.VxcUP)
VxcDOWN2 = collect(dis2.ksham.VxcDOWN)

pairs2 = filter(((i, j),) -> i <= n && j <= n, [(1, 1), (2, 2), (1, 2), (2, 3), (1, 4)])

@testset "GGA VxcUP/VxcDOWN vs finite difference of compute_exc_energy (nspin=2, PBE)" begin
    for (spin, Vxc) in ((1, VxcUP2), (2, VxcDOWN2))
        for (i, j) in pairs2
            D2p = copy(D2)
            D2m = copy(D2)
            D2p[i, j, spin] += h; i != j && (D2p[j, i, spin] += h)
            D2m[i, j, spin] -= h; i != j && (D2m[j, i, spin] -= h)
            Ep = AtomicKohnSham.compute_exc_energy(dis2, model2, D2p)
            Em = AtomicKohnSham.compute_exc_energy(dis2, model2, D2m)
            fd = (Ep - Em) / (2h)
            expected = i == j ? Vxc[i, i] : 2 * Vxc[i, j]
            @test isapprox(expected, fd; atol = 1e-6, rtol = 1e-4)
        end
    end
end
