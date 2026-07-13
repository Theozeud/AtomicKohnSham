# ===================================================================
#                          ALIGNED KEY = VALUE BLOCKS
# ===================================================================
"""
    _write_kv_block(io, indent, rows)

Print each `(label, value)` pair in `rows` as its own `label = value` line,
indented by `indent` spaces, with the `=` signs aligned to the longest label
in `rows`. Used throughout [`write_report`](@ref)/[`write_log_header`](@ref)
instead of hand-counted padding, which is easy to get wrong (and stays wrong
silently) whenever a label changes.
"""
function _write_kv_block(io::IO, indent::Int, rows)
    width = maximum(((label, _),) -> length(label), rows)
    pad = " "^indent
    for (label, value) in rows
        println(io, pad, rpad(label, width), " = ", value)
    end
end

# ===================================================================
#                     FUNCTIONAL / ALGORITHM NAMING
# ===================================================================
"""
Human-readable name of an exchange or correlation functional, used when
writing reports and logs.
"""
_functional_name(f) = hasproperty(f, :name) ? f.name : string(nameof(typeof(f)))
_functional_name(::NoFunctional) = "None"

# ===================================================================
#                       AUFBAU PARAMETER DUMP
# ===================================================================
function _write_aufbau(io::IO, aufbau::OptimizedAufbau)
    println(io, "    Aufbau = OptimizedAufbau")
    _write_kv_block(io, 6, [
        ("max_degen", aufbau.max_degen),
        ("tol", aufbau.tol),
        ("handle_degen_spin_1iter", aufbau.handle_degen_spin_1iter),
    ])
end

function _write_aufbau(io::IO, aufbau::SmearedAufbau)
    println(io, "    Aufbau = SmearedAufbau")
    _write_kv_block(io, 6, [("Temp", aufbau.Temp)])
end

function _write_aufbau(io::IO, aufbau::FrozenAufbau)
    println(io, "    Aufbau = FrozenAufbau")
    _write_kv_block(io, 6, [("n", aufbau.n)])
end

_write_aufbau(io::IO, aufbau::Aufbau) = println(io, "    Aufbau = $(nameof(typeof(aufbau)))")

# ===================================================================
#                     ALGORITHM PARAMETER DUMP
# ===================================================================
function _write_algorithm(io::IO, alg::ODA)
    println(io, "  Algorithm = ODA")
    _write_kv_block(io, 4, [
        ("scftol", alg.scftol),
        ("tinit", alg.tinit),
        ("frozen_t", alg.frozen_t),
        ("maxiter_ls", alg.maxiter_ls),
        ("abstol_ls", alg.abstol_ls),
        ("reltol_ls", alg.reltol_ls),
    ])
    _write_aufbau(io, alg.aufbau)
end

function _write_algorithm(io::IO, alg::SCFAlgorithm)
    println(io, "  Algorithm = $(nameof(typeof(alg)))")
    _write_kv_block(io, 4,
        [(string(name), getfield(alg, name)) for name in fieldnames(typeof(alg))])
end

# ===================================================================
#                      PARAMETERS BLOCK (shared)
# ===================================================================
"""
    _write_parameters(io, model, alg, lh, nh, basis, fem_integration_method)

Write the block describing the physical model, discretization, and SCF
algorithm of a computation. Shared between [`write_report`](@ref) (final
`results.txt`) and [`write_log_header`](@ref) (`log.txt` header), so that
both describe exactly the same run.
"""
function _write_parameters(io::IO, model::KSEModel, alg::SCFAlgorithm, lh::Int, nh::Int,
                          basis::FEMBasis, fem_integration_method::FEMIntegrationMethod)
    println(io, "PARAMETERS")
    println(io, "----------")
    println(io, "  Model")
    _write_kv_block(io, 4, [
        ("Z (nuclear charge)", model.Z),
        ("N (electrons)", model.N),
        ("Hartree coefficient", model.hartree),
        ("Exchange functional", _functional_name(model.exchange)),
        ("Correlation functional", _functional_name(model.correlation)),
        ("Spin channels", model.nspin),
    ])
    println(io)

    basis_value = "$(nameof(typeof(basis.generators))) (ordermax = $(basis.generators.ordermax))"
    mesh_value = basis.mesh.name * (isempty(basis.mesh.params) ? "" : " $(basis.mesh.params)")
    integration_value = string(nameof(typeof(fem_integration_method))) *
        (hasproperty(fem_integration_method, :npoints) ?
            " (npoints = $(fem_integration_method.npoints))" : "")

    println(io, "  Discretization")
    _write_kv_block(io, 4, [
        ("Angular momentum cutoff (lh)", lh),
        ("Orbitals per channel (nh)", nh),
        ("Basis", basis_value),
        ("Number of basis functions", length(basis)),
        ("Mesh", mesh_value),
        ("Rmax", last(basis.mesh)),
        ("Number of mesh points", length(basis.mesh)),
        ("Integration method", integration_value),
    ])
    println(io)

    _write_algorithm(io, alg)
end

# ===================================================================
#                        RESULTS.TXT REPORT
# ===================================================================
"""
    write_report(sol::KSESolution, path::AbstractString)

Write a self-contained text report of a converged (or stopped) Kohn–Sham
computation to `path`.

The report describes:
- the run's identity and outcome (name, success flag, number of iterations,
  final stopping criterion),
- every parameter used to define the computation (model, discretization,
  algorithm), read from `sol.context`,
- the final energy decomposition,
- the occupied orbitals with their energies and occupation numbers.

# Arguments
- `sol::KSESolution`: Solution returned by [`groundstate`](@ref).
- `path::AbstractString`: Destination file path (typically `"results.txt"`).
"""
function write_report(sol::KSESolution{T}, path::AbstractString) where T
    @unpack model, alg, lh, nh, basis, fem_integration_method = sol.context
    open(path, "w") do io
        _write_kv_block(io, 0, [
            ("Name", sol.name),
            ("Success", sol.success),
            ("niter", sol.niter),
        ])
        println(io, "Final stopping criterion = $(sol.stopping_criteria)")
        println(io, "Data type = $(T)")
        println(io)

        _write_parameters(io, model, alg, lh, nh, basis, fem_integration_method)
        println(io)

        println(io, "ENERGIES")
        println(io, "--------")
        _write_kv_block(io, 2, [
            ("Ekin (kinetic)", sol.energies.Ekin),
            ("Ecou (nuclear attraction)", sol.energies.Ecou),
            ("Ehar (Hartree)", sol.energies.Ehar),
            ("Eexc (exchange-correlation)", sol.energies.Eexc),
            ("Etot (total)", sol.energies.Etot),
        ])
        println(io)

        println(io, "OCCUPIED ORBITALS")
        println(io, "-----------------")
        _print_occupied(io, sol; styled = false)
        println(io)

        println(io, "SANITY CHECKS")
        println(io, "-------------")
        _print_sanity(io, sol; styled = false)
    end
    nothing
end
