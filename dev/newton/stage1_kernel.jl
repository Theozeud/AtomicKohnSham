# Stage 1: exchange-correlation kernels `eval_v2rho2` of the built-in functionals.
# Checked against central finite differences of `vrho` in Double64 (so that the
# finite differences are not limited by Float64 roundoff) and against Libxc
# `v2rho2`. The polarized PW92 kernel differs from Libxc by ~1e-7 relative, the
# same level as the known difference between the built-in PW92 `vrho` and Libxc.
# Run: julia stage1_kernel.jl
using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "..", "examples"))
using AtomicKohnSham, Libxc, DoubleFloats, Printf
using AtomicKohnSham: eval_v2rho2, eval_vrho, eval_vrho_up, eval_vrho_down

const T = Double64
const ρs = T[1e-6, 1e-4, 1e-2, 0.1, 1, 10, 1e3]
const ζs = T[-0.9, -0.3, 0, 0.2, 0.7, 0.95]
const h = T(1e-12)

for (name, lib) in (("Slater", :lda_x), ("PW92", :lda_c_pw))
    f1, f2 = BuiltinFunctional(lib; nspin = 1), BuiltinFunctional(lib; nspin = 2)
    x1, x2 = Functional(lib; n_spin = 1), Functional(lib; n_spin = 2)
    efd1, elx1, efd2, elx2 = 0.0, 0.0, 0.0, 0.0
    for ρ in ρs
        k = eval_v2rho2(f1, ρ)
        fd = (eval_vrho(f1, ρ + h*ρ) - eval_vrho(f1, ρ - h*ρ)) / (2h*ρ)
        lx = Libxc.evaluate(x1; rho = [Float64(ρ)], derivatives = 2:2).v2rho2[1]
        efd1 = max(efd1, Float64(abs(k - fd) / abs(fd)))
        elx1 = max(elx1, abs(Float64(k) - lx) / abs(lx))
        for ζ in ζs
            ρu, ρd = ρ * (1 + ζ) / 2, ρ * (1 - ζ) / 2
            k2 = eval_v2rho2(f2, ρu, ρd)
            hu, hd = h * ρu, h * ρd
            fd2 = ((eval_vrho_up(f2, ρu + hu, ρd) - eval_vrho_up(f2, ρu - hu, ρd)) / 2hu,
                   (eval_vrho_up(f2, ρu, ρd + hd) - eval_vrho_up(f2, ρu, ρd - hd)) / 2hd,
                   (eval_vrho_down(f2, ρu, ρd + hd) - eval_vrho_down(f2, ρu, ρd - hd)) / 2hd)
            lx2 = Libxc.evaluate(x2; rho = reshape(Float64[ρu, ρd], 2, 1),
                                 derivatives = 2:2).v2rho2[:, 1]
            sc = max(abs(fd2[1]), abs(fd2[3]))
            efd2 = max(efd2, Float64(maximum(abs(k2[c] - fd2[c]) for c in 1:3) / sc))
            elx2 = max(elx2, maximum(abs(Float64(k2[c]) - lx2[c]) for c in 1:3) / Float64(sc))
        end
    end
    @printf("%-7s unpolarized: vs FD %.1e, vs Libxc %.1e | polarized: vs FD %.1e, vs Libxc %.1e\n",
            name, efd1, elx1, efd2, elx2)
end
