"""
    line_search_energy(energy_kin0::T, energy_kin1::T,
                       energy_cou0::T, energy_cou1::T,
                       energy_har0::T, energy_har1::T,
                       energy_har01::T, energy_har10::T,
                       F::FT = zero(T);
                       nl::Bool = true,
                       maxiter::Int = 100,
                       abstol::T = eps(T)*1e4,
                       reltol::T = eps(T)*1e4) where {FT,T}

Compute the minimizer `t*` and the minimum value of the 1D energy along a mixing line
`t ∈ [0,1]`, for an energy model of the form

`E(t) = a t^2 + b t + c + F(t)`,

where the quadratic part comes from kinetic + Coulomb + Hartree contributions of two
states (e.g. density matrices / density operators) and `F(t)` is an optional nonlinear
correction (typically exchange–correlation, or any user-provided nonlinear term).

The coefficients are built from the inputs as:
- `A0 = energy_kin0 + energy_cou0`, `A1 = energy_kin1 + energy_cou1`
- `a = energy_har0 + energy_har1 - (energy_har01 + energy_har10)`
- `b = (A0 - A1) + (energy_har01 + energy_har10 - 2*energy_har1)`
- `c = A1 + energy_har1`

If `nl == true`, the function minimizes `E(t)` on `[0,1]` using Brent's method
(Optim.jl), with tolerances `abstol`, `reltol` and iteration cap `maxiter`.
If `nl == false`, it minimizes the quadratic part analytically on `[0,1]`.

This routine is used in ODA line-search steps, and can also be reused to optimize
fractional occupations of the last shell in case of degeneracy by choosing an
appropriate nonlinear term `F`.

# Arguments
- `energy_kin0, energy_kin1`: kinetic energies of the two endpoints.
- `energy_cou0, energy_cou1`: (external/ionic) Coulomb energies of the two endpoints.
- `energy_har0, energy_har1`: Hartree self-energies of the two endpoints.
- `energy_har01, energy_har10`: Hartree cross terms between endpoints.
- `F(t)`: callable nonlinear term (defaults to `zero(T)`).

# Returns
`(tmin, Emin)` where `tmin ∈ [0,1]` minimizes the objective and `Emin = E(tmin)`.
"""
function line_search_energy(energy_kin0::T, energy_kin1::T, energy_cou0::T, energy_cou1::T,
                           energy_har0::T, energy_har1::T, energy_har01::T, F::FT = zero(T);
                           nl::Bool = !(F isa Real), maxiter::Int = 100,
                           abstol::T = eps(T)*10^4, reltol::T = eps(T)*10^4) where {FT,T}
    A0 = energy_kin0 + energy_cou0
    A1 = energy_kin1 + energy_cou1
    H0 = energy_har0
    H1 = energy_har1
    H01 = energy_har01
    a = H0 + H1 - 2*H01
    b = (A0 - A1) + 2*(H01 - H1)
    c = A1 + H1
    @show nl
    if nl
        # Perform the optimization through the Brent method
        function f(t::DT) where DT
            a*t^2 + b*t + c + F(t)
        end
        X = LinRange(0,1,50)
        fX = f.(X)
        println(fX)

        res = optimize(f, zero(T), one(T), Brent();
                abs_tol=abstol, rel_tol=reltol, iterations=maxiter)
        fmin = res.minimum
        tmin = res.minimizer
        #@show (tmin, fmin)
        #=
        f0 = f(zero(T))
        f1 = f(one(T))
        return tmin, fmin
        if f0 <= fmin && f0 <= f1
            return zero(T), f0
        elseif f1 <= fmin && f1 <= f0
            return one(T), f1
        end
        =#
        return tmin, fmin
    else
        # No nonlinear term : optimization of a quadratic functional on [0,1]
        tmin = if a > 0
            clamp(-b/(2a), zero(T), one(T))
        else
            f0 = c
            f1 = c + b + a
            f0 <= f1 ? zero(T) : one(T)
        end
        fmin = c + b*tmin + a*tmin^2
        return tmin, fmin
    end
end
