# Stage 6: ODA on iron in reduced Hartree-Fock (degenerate 5s/3d Fermi level), run
# without stopping criterion, tracking the true self-consistency residual
# ‖D_out[D] - D‖: from iteration ~50 the line search returns t = 0 at every
# iteration and the state is frozen at a residual of 1.6e-6 (QCFON: 1e-11).
# Run: julia stage6_oda_stagnation.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", "examples"))
using AtomicKohnSham, LinearAlgebra, Printf
using AtomicKohnSham: empty_channels, diagonalize_channels!, density_matrix!, assemble_hamiltonian!

"""Every 25 iterations: ODA step length, mixing parameter t and the true residual ‖D_out[D] - D‖."""
struct Track end
function AtomicKohnSham.callback!(::Track, solver::KSESolver)
    it = solver.niter
    (it % 25 == 0 || it < 3) || return
    dis, model, c = solver.discretization, solver.model, solver.algcache
    nb = dis.Nₕ
    D = reshape(copy(c.D), nb, nb, 1)
    assemble_hamiltonian!(dis, model, c.D)
    ch = diagonalize_channels!(empty_channels(dis), dis, dis.ksham.H)
    for cc in ch
        nk = c.n[cc.l + 1, :, cc.σ]; cc.occ = findall(>(0), nk); cc.nocc = nk[cc.occ]
    end
    Dout = zeros(nb, nb, 1); density_matrix!(Dout, ch)
    n3d = c.n[3, 1, 1]
    @printf("it %4d  step=%.1e  t=%.3e  TRUE residual=%.2e  n(3d)=%.10f\n", it,
            solver.stopping_criteria, c.t, norm(Dout .- D), n3d)
end

model = RHF(Z = 26, N = 26)
basis = P1IntLegendreBasis(expmesh(0, 120, 25; s = 1.2); ordermax = 8)
dis = KSEDiscretization(basis, model; lh = 2, nh = 8)
alg = ODA(tinit = 0.6, scftol = 0.0, aufbau = OptimizedAufbau(max_degen = 2, tol = 1e-5))
groundstate(model, dis, alg; maxiter = 400, callback = Track())
