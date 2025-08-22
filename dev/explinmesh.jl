using AtomicKohnSham
function explinrange(a::Real, b::Real, n::Int;
                     rswitch::Real, s::Real=1.0, nlin::Int = 0, T=Float64)
    rexp= AtomicKohnSham.exprange(a, rswitch, n; T=T, s=s)
    if nlin > 0
        rlin = collect(T.(LinRange(rswitch, b, nlin+1)))
        return vcat(rexp[1:end-1], rlin)
    else
        return rexp
    end
end

function explinmesh(a::Real, b::Real, n::Int;
                    rswitch::Real, s::Real=1.0, nlin::Int, T=Float64)
    R = explinrange(a, b, n; rswitch=rswitch, nlin=nlin, s=s, T=T)
    AtomicKohnSham.Mesh(R, "Exp-Lin Mesh", (rswitch=rswitch, nlin=nlin, s=s))
end
