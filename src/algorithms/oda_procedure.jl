# COMPUTE THE MINIMUM AND THE MINIMIZER OF
#                       t ∈ [0,1] ↦ E(tγ₀ + (1-t)γ₁)
function find_minima_oda(energy_kin0::Real, energy_kin1::Real,
        energy_cou0::Real, energy_cou1::Real,
        energy_har0::Real, energy_har1::Real,
        energy_har01::Real, energy_har10::Real,
        D0::AbstractArray{<:Real}, D1::AbstractArray{<:Real}, tmpD::AbstractArray{<:Real},
        model::KSEModel, discretization::KSEDiscretization; occup::Bool = false)
    T = eltype(discretization)

    if has_exchcorr(model)
        # DEFINE THE OBJECTIV FUNCTION
        function f(t)
            @. tmpD = t*D0 .+ (1-t) * D1
            (t * (energy_kin0 + energy_cou0) + (1-t) * (energy_kin1 + energy_cou1) +
             t^2*energy_har0 + (1-t)^2 * energy_har1 +
             t * (1-t) * (energy_har01 + energy_har10) +
             compute_exchangecorrelation_energy(discretization, model, tmpD))
        end
        if occup
            #@show f.(LinRange(0,1,100))
        end
        # PERFORM THE OPTIMISATION THROUGHT THE GOLDEN SECTION's METHOD

        res = optimize(f, zero(T), one(T), Brent();
                abs_tol=1e-10, rel_tol=1e-10, iterations=100)
        fmin = res.minimum
        tmin = res.minimizer
        f0 = f(zero(T))
        f1 = f(one(T))
        return tmin, fmin
        if f0 <= fmin && f0 <= f1
            return zero(T), f0
        elseif f1 <= fmin && f1 <= f0
            return one(T), f1
        end
        return tmin, fmin
    else
        # DEFINE THE OBJECTIV FUNCTION
        function f2(t)
            t * (energy_kin0 + energy_cou0) + (1-t) * (energy_kin1 + energy_cou1) +
            t^2*energy_har0 + (1-t)^2 * energy_har1 +
            t * (1-t) * (energy_har01 + energy_har10)
        end

        # PERFORM THE OPTIMISATION THROUGHT THE GOLDEN SECTION's METHOD
        res = optimize(f2, zero(T), one(T), Brent())

        return res.minimizer, res.minimum

    end
end
