# STRUCTURE TO CHANGE THE METHOD DURING THE PROCESS OF MINIMIZATION

#####################################################################
#                          STRUCTURE
#####################################################################

mutable struct CombinedMethod{typeMethod1 <: SCFAlgorithm, typeMethod2 <: SCFAlgorithm} <: SCFAlgorithm
    method1::typeMethod1
    method2::typeMethod2
    state::Int
    iter_switch::Int
end

mutable struct CombinedCache{typeCache1 <: SCFCache, typeCache2 <: SCFCache} <: SCFCache
    cache1::typeCache1
    cache2::typeCache2
end

function create_cache_alg(method::CombinedMethod, discretization::KSEDiscretization)
    cache1 = create_cache_alg(method.method1, discretization)
    cache2 = create_cache_alg(method.method2, discretization)
    CombinedCache{typeof(cache1), typeof(cache2)}(cache1, cache2)
end


#####################################################################
#                              STEPS
#####################################################################

function loopheader!(cache::CombinedCache, method::CombinedMethod)
    if solver.iter == method.iter_switch
        method.state = 2
        switch!(cache2, cache1)
    end
    if method.state == 1
        loopheader!(cache.cache1, method.method1)
    else
        loopheader!(cache.cache2, method.method2)
    end
end

function performstep!(cache::CombinedCache, method::CombinedMethod, solver::KSESolver) 
    if method.state == 1
        performstep!(cache.cache1, method.method1, solver)
    else
        performstep!(cache.cache2, method.method2, solver)
    end
end

function loopfooter!(cache::CombinedCache, method::CombinedMethod)
    if method.state == 1
        loopfooter!(cache.cache1, method.method1)
    else
        loopfooter!(cache.cache2, method.method2)
    end
end


function monitor(cache::CombinedCache, method::CombinedMethod)
    if method.state == 1
        monitor(cache.cache1, method.method1)
    else
        monitor(cache.cache2, method.method2)
    end
end


#####################################################################
#                              SWITCH
#####################################################################
#=
function switch!(cache2::CacheQuadratic, cache1::RCACache)


end
=#