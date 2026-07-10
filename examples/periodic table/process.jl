using AtomicKohnSham
using LinearAlgebra
using Base.Threads
using Libxc


function predic_nb_iter(init::Real, tol::Int)
    Int(floor(1 + (tol*log(10) + log(init))/log(2)))
end


function compute_maximum_ionization(z::Int, spin_sym::Bool)

    # INITIALIZATION OF PARAMETERS
    binf = z
    bsup = z + 1.2
    mid = (bsup + binf)/2
    precision = 2
    tol = 10.0^(-precision)
    lh = Int(min(ceil(z/20),3))
    nh = min(max(z,2),10)
    scftol = 1e-6
    Nmesh = 20

    # SELECTION OF RMAX
    problem = AtomProblem(; z = z,
                            N = z,
                            ex = Functional(:lda_x, n_spin = 2-spin_sym),
                            ec = Functional(:lda_c_pw, n_spin = 2-spin_sym),
                            lh = lh,
                            nh = nh,
                            Nmesh = Nmesh,
                            Rmax = 10000,
                            typemesh = expmesh,
                            optsmesh = (s = 1.3,),
                            typebasis = P1IntLegendreBasis,
                            optsbasis = (ordermax = 10,),
                            integration_method = GaussLegendre,
                            optsintegration = (npoints = 2000,),
                            alg = ODA(0.4),
                            scftol = scftol,
                            maxiter = 100,
                            degen_tol = 1e-2)

    sol = groundstate(problem)
    X = expmesh(0,10000,10*Nmesh;s=1.5).points
    ρX2 = eval_density(sol, X) .* X.^2
    i = findfirst(x->x<1e-20, ρX2)
    Rmax = X[i]*4

    # DICHOTOMIE TO FIND THE MAXIMUM OF IONIZATION
    while abs(bsup - binf) > tol
        mid = (bsup + binf)/2
        # CREATION OF THE PROBLEM
        lh = Int(min(ceil(mid/20),3))
        problem = AtomProblem(; z = z,
                                N = mid,
                                ex = Functional(:lda_x, n_spin = 2-spin_sym),
                                ec = Functional(:lda_c_pw, n_spin = 2-spin_sym),
                                lh = lh,
                                nh = nh,
                                Nmesh = Nmesh,
                                Rmax = Rmax,
                                typemesh = expmesh,
                                optsmesh = (s = 1.3,),
                                typebasis = P1IntLegendreBasis,
                                optsbasis = (ordermax = 10,),
                                integration_method = GaussLegendre,
                                optsintegration = (npoints = 1000,),
                                alg = ODA(0.4),
                                scftol = scftol,
                                maxiter = 50,
                                degen_tol = 1e-2)
        # RESOLUTION
        sol = groundstate(problem)

        # CHEK STABILITY
        ϵlast = last(sol.occupation_number)[2]
        if ϵlast > 0 || sol.stopping_criteria > scftol
            bsup = mid
        else
            binf = mid
        end

    end
    return round(mid;digits=precision)
end


function max_ionization_all_elements(zmin::Int=1, zmax::Int=118)
    result = zeros(zmax-zmin+1)
    Threads.@threads for z ∈ zmin:zmax
        @time "z = $z" begin
            ioni_nosym  = compute_maximum_ionization(z,false) - z
            ioni_sym    = compute_maximum_ionization(z,true) - z
        end
        result[z-zmin+1] = max(ioni_nosym,ioni_sym)
    end
    result
end
