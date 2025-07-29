#  CONVERGENCE WITH NMESH

function convergenceNmesh(vecNmesh::AbstractVector, problems; nums = [1])
    ΔΛ = Dict()
    ΔU = Dict()
    for problem in problems
        @unpack T, name = problem
        println(name)
        Δλ = zeros(T, length(vecNmesh), length(nums))
        @inbounds for i in eachindex(vecNmesh)
            newprob = HydrogenoidProblem(problem; Nmesh = vecNmesh[i], nλ = nums, nU = nums)
            @time "Nmesh = $(vecNmesh[i])" sol = eigvals_hydro(newprob)
            Δλ[i, :] .= sol.Δλ
        end
        ΔΛ[name] = Δλ
    end
    HydrogenoidConvergenceNmesh(problems, vecNmesh, ΔΛ, ΔU, nums)
end

#  CONVERGENCE WITH RMAX

function convergenceRmax(vecRmax::AbstractVector, problems; nums = [1])
    ΔΛ = Dict()
    ΔU = Dict()
    for problem in problems
        @unpack T, name = problem
        println(name)
        Δλ = zeros(T, length(vecNmesh), length(nums))
        @inbounds for i in eachindex(vecRmax)
            newprob = HydrogenoidProblem(problem; Rmax = vecRmax[i], nλ = nums, nU = nums)
            @time "Rmax = $(vecRmax[i])" sol = eigvals_hydro(newprob)
            Δλ[i, :] = sol.Δλ
        end
        ΔΛ[name] = Δλ
    end
    HydrogenoidConvergenceNmesh(problems, vecRmax, ΔΛ, ΔU, nums)
end
