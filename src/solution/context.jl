"""
    KSEContext

Immutable context describing a Kohn–Sham calculation.

Stores the physical model, discretization parameters, basis, and integration
method required to interpret and post-process a converged Kohn–Sham solution.
This structure is independent of the SCF solver state and contains no mutable
workspaces.
"""
struct KSEContext{M<: KSEModel, B<:FEMBasis, F<: FEMIntegrationMethod}
    model::M
    alg::A
    lh::Int
    nh::Int
    basis::B
    fem_integration_method::F
end
