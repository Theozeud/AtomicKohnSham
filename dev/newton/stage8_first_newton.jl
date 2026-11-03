# Stage 8: the first Newton step of QCFON, which does not decrease the residual
# (e.g. 8.6e-3 -> 1.3e-2 on chlorine near its threshold). From the ODA iterate D₀
# at the switch and the Newton step δD (read in a callback after the first Newton
# iteration: D₀ = cache.Dstart, δD = cache.Dn - D₀), the residual
# F(D) = D_out[D] - D (occupations of the current structure) is evaluated along
# D₀ + t δD. With an exact Jacobian and an exact linear solve,
# F(D₀ + t δD) = (1 - t) F(D₀) + O(t²): the deviation from the linear prediction
# scales as t² if the step is right and the first step fails by nonlinearity, and
# has a t term otherwise (wrong Jacobian or inexact GMRES).
# Run: julia stage8_first_newton.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", "examples"))
using AtomicKohnSham, LinearAlgebra, Printf
using AtomicKohnSham: diagonalize_channels!, density_matrix!, assemble_hamiltonian!,
    QCFONCache, _spin

"""F(D) = D_out[D] - D with the occupations of `channels`, per spin block."""
function residual(dis, model, channels, D)
    nspin = dis.nspin
    assemble_hamiltonian!(dis, model, nspin == 1 ? D[:, :, 1] : D)
    ch = diagonalize_channels!(deepcopy(channels), dis, dis.ksham.H)
    Dout = similar(D)
    density_matrix!(Dout, ch)
    Dout .- D, ch
end

"""Line study after the first Newton step; `done` prevents repeating it."""
mutable struct FirstStep
    done::Bool
end
function AtomicKohnSham.callback!(f::FirstStep, solver::KSESolver)
    c = solver.algcache
    (f.done || !(c isa QCFONCache) || !c.newton || iszero(c.prevcrit)) && return
    f.done = true
    dis, model = solver.discretization, solver.model
    nb, nspin = dis.Nₕ, dis.nspin
    D0 = reshape(copy(c.Dstart), nb, nb, nspin)
    δD = c.Dn .- D0
    F0, _ = residual(dis, model, c.channels, D0)
    @printf("  ‖F(D₀)‖ = %.3e   ‖δD‖ = %.3e   fractional orbitals: %d\n",
            norm(F0), norm(δD), length(c.partial))
    @printf("  %8s %12s %14s %10s\n", "t", "‖F(D₀+tδD)‖", "‖F-(1-t)F₀‖", "/t²")
    for t in (1e-3, 1e-2, 1e-1, 0.25, 0.5, 0.75, 1.0)
        Ft, _ = residual(dis, model, c.channels, D0 .+ t .* δD)
        dev = norm(Ft .- (1 - t) .* F0)
        @printf("  %8.3f %12.3e %14.3e %10.3e\n", t, norm(Ft), dev, dev / t^2)
    end
    nothing
end

function study(name, Z, N, nspin, Rmax, nmesh)
    ex = BuiltinFunctional(:lda_x; nspin = nspin)
    ec = BuiltinFunctional(:lda_c_pw; nspin = nspin)
    model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)
    basis = P1IntLegendreBasis(expmesh(0, Rmax, nmesh; s = 1.2); ordermax = 10)
    dis = KSEDiscretization(basis, model; lh = 2, nh = 4)
    println(name)
    groundstate(model, dis, QCFON(scftol = 1e-11); maxiter = 200, callback = FirstStep(false))
end

study("Ne, LDA", 10, 10, 1, 100, 25)
study("Cl, LSDA", 17, 17, 2, 100, 25)
study("Cl, LSDA, N = 17.9864", 17, 17.9864, 2, 200, 30)
