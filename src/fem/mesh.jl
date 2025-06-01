#-------------------------------------
#            MESH STRUCTURE
#-------------------------------------
struct Mesh{T <: Real,
            S <: AbstractString,
            P <: AbstractVector,
            T <: NamedTuple} 
    name::S
    points::P
    params::T
    function Mesh(points::AbstractVector{T}, name::AbstractString = "", params::NamedTuple = NamedTuple{}())  where T <: Real
        new{T, typeof(name), typeof(points), typeof(params)}(points, name, params)
    end
end

#-------------------------------------
#                   API
#-------------------------------------

function Base.show(io::IO, m::Mesh)
    println(io, "Mesh: \"", m.name, "\"")
    if !isempty(m.params)
        println(io, "Parameters:")
        for (k,v) in pairs(m.params)
            println(io, "  $k => $v")
        end
    end
    println(io, "  Number of points: ", length(m.points))
    println(io, "  First point: ", first(m.points))
    println(io, "  Last point: ", last(m.points))
    println(io, "  Points: ", m.points)
end

@inline Base.eltype(::Mesh{T}) where T = T
@inline Base.eachindex(m::Mesh) = eachindex(m.points)
@inline Base.firstindex(m::Mesh) = firstindex(m.points)
@inline Base.lastindex(m::Mesh) = lastindex(m.points)
@inline Base.first(m::Mesh) = first(m.points)
@inline Base.last(m::Mesh) = last(m.points)
@inline Base.getindex(m::Mesh, n::Int) = m.points[n]
@inline Base.getindex(m::Mesh, ur::UnitRange{Int64}) = m.points[ur]
@inline Base.length(m::Mesh) = length(m.points)
@inline Base.size(m::Mesh) = size(m.points)
@inline Base.iterate(m::Mesh, state = 1) = state > length(m) ? nothing : (m[state],state+1)
@inline cellrange(m::Mesh) = firstindex(m):lastindex(m)-1


@inline function findindex(m::Mesh, x::Real)
    if x ≤ m[end]
        return searchsortedlast(m.points, x)
    else
        return lastindex(m)+1
    end
end


#-------------------------------------
#           LINEAR MESH
#-------------------------------------

linmesh(a::Real, b::Real, n::Int; T::Type = Float64) = Mesh(T.(LinRange(a,b,n)), "Linear Mesh")


#-------------------------------------
#           GEOMETRIC MESH
#-------------------------------------

function geometricrange(a::Real ,b::Real ,n::Int; T::Type = Float64, s::Real)
    R = zeros(T,n)
    R[1] = a
    hn = (one(T)-T(s))/(one(T) - T(s)^(n-1))*(b-a)
    H = zeros(T,n-1)
    H[end] = hn
    for i ∈ n-1:-1:2
        H[i-1]= T(s) * H[i]
    end 
    for i ∈ 2:n
        R[i] = R[i-1] + H[i-1]
    end
    R
end

geometricmesh(a::Real, b::Real ,n::Int; T = Float64, s::Real) = Mesh(geometricrange(a,b,n; T = T, s = s), "Geometric Mesh", (s = s,))


#-------------------------------------
#           POLYNOMIAL MESH
#-------------------------------------

function polynomialrange(a::Real, b::Real, n::Int; T::Type = Float64, s::Real)
    R = zeros(T,n)
    R[1] = T(a)
    R[end] = T(b)
    for i ∈ 2:n-1
        R[i] = ((i-1)/(n-1))^s * (T(b)-T(a)) + T(a)  
    end
    R
end

polynomialmesh(a::Real, b::Real ,n::Int; T = Float64, s) = Mesh(polynomialrange(a,b,n; T = T, s = s), "Polynomial Mesh", (s = s,))


#-------------------------------------
#          EXPONENTIAL MESH
#-------------------------------------

function exprange(a::Real, b::Real, n::Int; T::Type = Float64, s::Real)
    R = zeros(T,n)
    R[1] = T(a)
    R[end] = T(b)
    for i ∈ 2:n-1
        pow = ((i-1)/(n-1))^s
        R[i] = (1 + b-a)^pow - 1 + a
    end
    R
end

expmesh(a::Real, b::Real ,n::Int; T = Float64, s) = Mesh(exprange(a,b,n; T = T, s = s), "Exponential Mesh", (s = s,))