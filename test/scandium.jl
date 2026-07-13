function run_solve(ex::TEX, ec::TEC,
                   Rmax::Real, Nmesh::Int, ordermax::Int,
                   maxiter::Int, scftol::Real) where {TEX, TEC}
    # Creation of the model
    model = KSEModel(; Z = 21, N = 21, ex = ex, ec = ec)

    # Creation of the mesh
    mesh = expmesh(0, Rmax, Nmesh; s = 1.2)

    # Creation of the Basis
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)

    # Creation of the discretization
    dis = KSEDiscretization(basis, model; lh = 2, nh=5)

    # Creation of the algorithm
    aufbau = OptimizedAufbau(max_degen = 2)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = scftol)

    # Resolution of the problem
    sol = groundstate(model, dis, alg; maxiter = maxiter)
    return sol
end

# See hydrogen.jl for the documented version of this helper.
function check_sanity(sol; tol::Real = 1e-6, virial_tol::Union{Real, Nothing} = nothing)
    sc = sol.sanity
    @test abs(sc.energy_balance) < tol
    @test abs(sc.electron_count_error) < tol
    @test abs(sc.density_norm_error) < tol
    @test abs(sc.orthonormality_error) < tol
    if !isnothing(virial_tol)
        @test abs(sc.virial_ratio - 1) < virial_tol
    end
end


@testset "Scandium RHF Spin Nopoloarized" begin

    # MODEL
    ex = NoFunctional(n_spin=1)
    ec = NoFunctional(n_spin=1)
    # DISCRETIZATION : P1IntLegendreBasis + expmesh
    Rmax        = 2000
    Nmesh       = 40
    ordermax    = 10
    # ALGORITHM : ODA + Optimized Aufbau
    maxiter     = 100
    scftol      = 1e-10



    @test_nowarn run_solve(ex, ec, Rmax, Nmesh, ordermax, maxiter, scftol)
    sol = run_solve(ex, ec, Rmax, Nmesh, ordermax, maxiter, scftol)

    @show sol
    @show sol.energies.Ekin

    @test abs(sol.energies.Etot + 722.5838224454205) < 1e-9

    occ1S = sol.occupied[1]
    @test occ1S[1] == "1s"
    @test abs(occ1S[2] + 154.35865) < 1e-5
    @test occ1S[3] == 2.0

    occ2S = sol.occupied[2]
    @test occ2S[1] == "2s"
    @test abs(occ2S[2] + 15.78538) < 1e-5
    @test occ2S[3] == 2.0

    occ2p = sol.occupied[3]
    @test occ2p[1] == "2p"
    @test abs(occ2p[2] + 12.74151) < 1e-5
    @test occ2p[3] == 6.0

    occ3S = sol.occupied[4]
    @test occ3S[1] == "3s"
    @test abs(occ3S[2] + 1.69002) < 1e-5
    @test occ3S[3] == 2.0

    occ3p = sol.occupied[5]
    @test occ3p[1] == "3p"
    @test abs(occ3p[2] + 0.96964) < 1e-5
    @test occ3p[3] == 6.0

    occ4S = sol.occupied[6]
    @test occ4S[1] == "4s"
    @test abs(occ4S[2] + 0.08646) < 1e-5
    @test occ4S[3] == 2.0

    occ4p = sol.occupied[7]
    @test occ4p[1] == "4p"
    @test abs(occ4p[2] + 0.00262) < 1e-5

    occ3d = sol.occupied[8]
    @test occ3d[1] == "3d"
    @test abs(occ3d[2] + 0.00262) < 1e-5
    @test abs(occ3d[3]/5 - 0.0056) < 1e-4

    # No exchange-correlation here (only nuclear Coulomb + Hartree, both ~1/r):
    # the virial theorem Ekin = -Etot applies, up to the Rmax-truncation defect.
    check_sanity(sol; tol = 1e-5, virial_tol = 1e-5)
end
