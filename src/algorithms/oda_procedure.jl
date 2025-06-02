# COMPUTE THE MINIMUM AND THE MINIMIZER OF 
#                       t ∈ [0,1] ↦ E(tγ₀ + (1-t)γ₁)
function find_minima_oda(energy_kin0::Real, energy_kin1::Real, 
                         energy_cou0::Real, energy_cou1::Real, 
                         energy_har0::Real, energy_har1::Real, 
                         energy_har01::Real, energy_har10::Real,
                         D0::AbstractArray{<:Real}, D1::AbstractArray{<:Real}, tmpD::AbstractArray{<:Real}, 
                         model::KSEModel, discretization::KSEDiscretization)
    
    T = eltype(discretization)

    if isthereExchangeCorrelation(model)
        # DEFINE THE OBJECTIV FUNCTION
        function f(t)        
            @. tmpD = t*D0 .+ (1-t) * D1
            (t * (energy_kin0 + energy_cou0) + (1-t) * (energy_kin1 + energy_cou1) + 
            t^2*energy_har0 + (1-t)^2 * energy_har1 + t*(1-t) * (energy_har01 + energy_har10) +
            compute_exchangecorrelation_energy(discretization, model, tmpD))
        end

        # PERFORM THE OPTIMISATION THROUGHT THE GOLDEN SECTION's METHOD
        res = optimize(f, zero(T), one(T), GoldenSection())

        return res.minimizer, res.minimum
    else
        # DEFINE THE OBJECTIV FUNCTION
        function f2(t)        
            t * (energy_kin0 + energy_cou0) + (1-t) * (energy_kin1 + energy_cou1) + 
            t^2*energy_har0 + (1-t)^2 * energy_har1 + t*(1-t) * (energy_har01 + energy_har10)
        end

        # PERFORM THE OPTIMISATION THROUGHT THE GOLDEN SECTION's METHOD
        res = optimize(f2, zero(T), one(T), GoldenSection())

        return res.minimizer, res.minimum
    end
end