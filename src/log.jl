#--------------------------------------------------------------------
#                               LogConfig
#--------------------------------------------------------------------
struct LogConfig{C}
    stopping_criteria::Bool
    energies::Bool
    methodlogconfig::C
    function LogConfig(;stopping_criteria = true, energies = true, kwargs...)
        new{typeof(kwargs)}(stopping_criteria, energies, kwargs)
    end
end


function Base.show(io::IO, lc::LogConfig)
    println(io, "LogConfig:")
    println(io, "  stopping_criteria      = ", lc.stopping_criteria)
    println(io, "  energies               = ", lc.energies)
    for k ∈ keys(lc.methodlogconfig)
        mk = lc.methodlogconfig[k]
        smk = string(mk)
        println(io, "  ",smk,repeat(" ",23 - length(mk)),"= ", mk)
    end
end


#--------------------------------------------------------------------
#                               LogBook
#--------------------------------------------------------------------
abstract type AbstractLogBook end 


struct LogBook{T<:Real, L<:AbstractLogBook} <: AbstractLogBook 
    config::LogConfig
    stopping_criteria::Vector{T}
    energies::Dict{Symbol, Vector{T}}
    methodlog::L

    function LogBook(config, T, energies::Dict, alg::SCFAlgorithm)
        stopping_criteria   = T[]               
        energies            = Dict(k => T[] for k in keys(energies))                
        methodlog           = create_logbook(alg)                 
        new{T, typeof(methodlog)}(config, stopping_criteria, energies, methodlog)
    end
end


function Base.getproperty(logbook::LogBook, s::Symbol)
    if s ∈ fieldnames(LogBook)
        getfield(logbook, s)
    elseif s ∈ propertynames(getfield(logbook, :methodlog))
        getfield(getfield(logbook, :methodlog), s)
    else
        throw(ErrorException("type LogBook has no field $(s)"))
    end
end
