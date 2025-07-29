struct SlaterXα <: BuiltinFunctional
    identifier::Symbol
    n_spin::Int
    name::String
    family::Symbol
    function SlaterXα(n_spin::Int)
        new(:lda_x, n_spin, "Slater Exchange", :lda)
    end
end

eval_zk(::SlaterXα, ρ::Real) = - 3/4 * (3/π)^(1/3) * ρ^(1/3)
eval_vrho(::SlaterXα, ρ::Real) = - (3/π)^(1/3) * ρ^(1/3)
