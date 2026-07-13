# Regenerates test/reference_data/fem_matrices.jls, the reference values the
# "FEM Matrices" testset in fem.jl compares fresh computations against.
#
# Not run as part of the test suite. Run manually after a deliberate change to
# the basis/integration code that should change these matrices/tensors:
#
#     julia --project=. test/generate_fem_reference_data.jl

using AtomicKohnSham
using Serialization

const REFERENCE_FILE = joinpath(@__DIR__, "reference_data", "fem_matrices.jls")

# Must match the basis construction in the "FEM Matrices" testset in fem.jl exactly.
function reference_basis()
    a, b, n, T, s = 0.0, 1.0, 3, Float64, 0.9
    mesh = polynomialmesh(a, b, n; T = T, s = s)
    P1IntLegendreBasis(mesh, T; ordermax = 2)
end

function compute_reference_data()
    basis = reference_basis()
    Dict(
        "M0" => mass_matrix(basis),           # overlap matrix
        "M1" => mass_matrix(basis, -1),       # overlap matrix, weight 1/x
        "M2" => mass_matrix(basis, -2),       # overlap matrix, weight 1/x^2
        "A" => stiffness_matrix(basis),       # stiffness matrix
        "T0" => mass_tensor(basis),           # mass tensor
        "T1" => mass_tensor(basis, -1),       # mass tensor, weight 1/x
        "T2" => mass_tensor(basis, -2),       # mass tensor, weight 1/x^2
    )
end

mkpath(dirname(REFERENCE_FILE))
serialize(REFERENCE_FILE, compute_reference_data())
println("Wrote reference data to $REFERENCE_FILE")
