abstract type SCFAlgorithm end
abstract type SCFCache end
abstract type SCFSolution end

#=
loopheader!(::SCFCache, ::SCFMethod, ::KohnShamSolver)      = nothing
performstep!(::SCFCache, m::SCFMethod, ::KohnShamSolver)    = @warn "No performstep for the method $(typeof(m))"
loopfooter!(::SCFCache, ::SCFMethod, ::KohnShamSolver)      = nothing
monitor(::SCFCache, ::SCFMethod, ::KohnShamSolver)          = nothing
register!(::SCFCache, ::SCFMethod, ::KohnShamSolver)        = nothing
create_cache_method(m::SCFMethod, 
                        ::KohnShamDiscretization)               = @warn "No creation of cache for the method $(typeof(m))"
switch!(c2::SCFCache, c1::SCFCache)                         = @warn "No way to init $(typeof(c2)) from of cache for the method ($(typeof(c1))"
=#
