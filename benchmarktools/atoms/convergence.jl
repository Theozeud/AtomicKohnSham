#--------------------------------------------------------------------
#                   ANALYSING CONVERGENCE IN NMESH
#--------------------------------------------------------------------
struct AtomConvergenceNmesh
    vecNmesh::Any        # Set of Nmesh used
    Error::Any           # Dict of errors on orbitals energy : for each problem,
    # there is a vector of errors depending on Nmesh
    num::Any             # Number of orbitals energy used 
end

function convergenceNmesh(vecNmesh::AbstractVector, problems; nums = [1])
    Errors = Dict()
    for problem in problems
        @unpack T, name = problem
        println(name)
        ϵ = zeros(T, length(vecNmesh), length(nums))
        @inbounds for i in eachindex(vecNmesh)
            newprob = AtomProblem(problem; Nmesh = vecNmesh[i])
            @time "Nmesh = $(vecNmesh[i])" sol = groundstate(newprob)
            ϵ[i, :] .= sol.orbitals_energy[nums]
        end
        error = zeros(T, length(vecNmesh)-1, length(nums))
        for i in axes(error, 1)
            error[i, :] = abs.(ϵ[i, :] .- ϵ[end, :])
        end
        Errors[name] = error
    end
    AtomConvergenceNmesh(vecNmesh, Errors, nums)
end

function convergencePlotNmesh(sols::AtomConvergenceNmesh, nums = firstindex(sols.num))
    plt = plot(
        size = (650, 500), margin = 0.5Plots.cm, legend = :outertopright, yaxis = :log,
        legendfontsize = 12,
        titlefontsize = 12,
        guidefontsize = 12,
        tickfontsize = 12)
    xlabel!(plt, "Nmesh")
    ylabel!(plt, "Error")
    title!(plt, "Convergence Plot")
    for key in keys(sols.Error)
        for num in nums
            plot!(plt, sols.vecNmesh[1:(end - 1)], sols.Error[key][:, num],
                lw = 4, label = key*"-$num", markershape = :x, markersize = 10)
        end
    end
    plt
end

#--------------------------------------------------------------------
#                   ANALYSING CONVERGENCE IN RMAX
#--------------------------------------------------------------------
struct AtomConvergenceRmax
    vecRmax::Any         # Set of Nmesh used
    Error::Any           # Dict of errors on orbitals energy : for each problem,
    # there is a vector of errors depending on Rmax
    num::Any             # Number of orbitals energy used 
end
function convergenceRmax(vecRmax::AbstractVector, problems; nums = [1])
    Errors = Dict()
    for problem in problems
        @unpack T, name = problem
        println(name)
        ϵ = zeros(T, length(vecRmax), length(nums))
        @inbounds for i in eachindex(vecRmax)
            newprob = AtomProblem(problem; Rmax = vecRmax[i])
            @time "Rmax = $(vecRmax[i])" sol = groundstate(newprob)
            ϵ[i, :] .= sol.orbitals_energy[nums]
        end
        error = zeros(T, length(vecRmax)-1, length(nums))
        for i in axes(error, 1)
            error[i, :] = abs.(ϵ[i, :] .- ϵ[end, :])
        end
        Errors[name] = error
    end
    AtomConvergenceRmax(vecRmax, Errors, nums)
end

function convergencePlotRmax(sols::AtomConvergenceRmax, nums = first(sols.num))
    plt = plot(
        size = (650, 500), margin = 0.5Plots.cm, legend = :outertopright, yaxis = :log,
        legendfontsize = 12,
        titlefontsize = 12,
        guidefontsize = 12,
        tickfontsize = 12)
    xlabel!(plt, "Rmax")
    ylabel!(plt, "Error")
    title!(plt, "Convergence Plot")
    for key in keys(sols.Error)
        for num in nums
            plot!(plt, sols.vecRmax[1:(end - 1)], sols.Error[key][:, num],
                lw = 4, label = key*"-$num", markershape = :x, markersize = 10)
        end
    end
    plt
end
