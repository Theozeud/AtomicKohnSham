#cd("./examples")
using Pkg
Pkg.activate(".")
using AtomicKohnSham
using Libxc

# ===================== STEP 0 : LOG =====================#
VERBOSE = 3
WRITE_LOG = true
LOGFILE = joinpath(@__DIR__, "scf.log")

# ===================== STEP 1 : PARAMETERS =====================#
# MODEL
Z = 3
N = 3
ex = Functional(:lda_x; n_spin=2)
ec = Functional(:lda_c_pw, n_spin=2)

# DISCRETIZATION : P1IntLegendreBasis + expmesh
Rmax        = 500
Nmesh       = 60
ordermax    = 10
lh          = 2

# ALGORITHM : ODA + Optimized Aufbau
maxiter     = 100
scftol      = 1e-11
ndict = Dict(
    "1sUP" => 1,
    "1sDOWN" => 1,
    "2sUP" => 0.975,
    "2sDOWN" => 0.025
)

# PLOT
Rmax_plot   = 100
N_plot      = 10000
lplot       = lh

# ===================== STEP 2 : SOLVE =====================#
function run_solve(Z::Real, N::Real, ex::TEX, ec::TEC,
                   Rmax::Real, Nmesh::Int, ordermax::Int, lh::Int,
                   maxiter::Int, ndict::Dict{String,<:Real}, scftol::Real) where {TEX, TEC}
    # Creation of the model
    model = KSEModel(; Z = Z, N = N, ex = ex, ec = ec)

    # Creation of the mesh
    mesh = expmesh(0, Rmax, Nmesh; s = 1.2)

    # Creation of the Basis
    basis = P1IntLegendreBasis(mesh; ordermax = ordermax)

    # Creation of the discretization
    dis = KSEDiscretization(basis, model; lh = lh, nh=10, fem_integration_method=GaussLegendre(basis, 2000))

    # Creation of the algorithm
    aufbau = FrozenAufbau(ndict)
    alg = ODA(tinit = 0.6, aufbau = aufbau, scftol = scftol)

    # Resolution of the problem
    sol = if WRITE_LOG
        open(LOGFILE, "w") do io
            redirect_stdout(io) do
                groundstate(model, dis, alg; maxiter = maxiter, verbose = VERBOSE)
            end
        end
    else
        sol = groundstate(model, dis, alg; maxiter = maxiter, verbose = VERBOSE)
    end
    return sol
end

sol = run_solve(Z, N,ex, ec, Rmax, Nmesh, ordermax, lh, maxiter, ndict, scftol)
