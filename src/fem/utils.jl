#####################################################################
#                          INTERSECTION VECTOR
#####################################################################

function find_intersection!(I::Vector{Int}, A::NTuple{N,Int}, B::NTuple{M,Int}) where {N,M}
    count = 0
    for a ∈ A
        if a ∈ B
            count += 1
            I[count] = a
        end
    end
    count
end


function find_intersection!(I::Vector{Int}, A::Int, B::NTuple{M,Int}) where {M}
    for b ∈ B
        if b == A
            I[1] = A
            return 1
        end
    end
    return 0
end

function find_intersection!(I::Vector{Int}, A::NTuple{N,Int}, B::Int) where {N}
    for a ∈ A
        if a == B
            I[1] = B
            return 1
        end
    end
    return 0
end

function find_intersection!(I::Vector{Int}, A::Int, B::Int) 
    if A == B
        I[1] = A
        return 1
    end
    return 0
end

#####################################################################
#                          INTERSECTION 2 ELEMENTS
#####################################################################

function find_intersection_indices!(I::Matrix{Int}, A::Int, B::Int)
    if A == B
        I[1,1] = 1
        I[2,1] = 1
        return 1
    end
    return 0
end


function find_intersection_indices!(I::Matrix{Int}, A::NTuple{N,Int}, B::Int) where {N}
    for i ∈ eachindex(A)
        if A[i] == B
            I[1,1] = i
            I[2,1] = 1
            return 1
        end
    end
    return 0
end

function find_intersection_indices!(I::Matrix{Int}, A::Int, B::NTuple{N,Int}) where {N}
    for i ∈ eachindex(B)
        if B[i] == A
            I[1,1] = 1
            I[2,1] = i
            return 1
        end
    end
    return 0
end

function find_intersection_indices!(I::Matrix{Int}, A::NTuple{N,Int}, B::NTuple{M,Int}) where {N,M}
    count = 0
    for a ∈ A
        if a ∈ B
            count += 1
            I[1,count] = findfirst(==(a),A)
            I[2,count] = findfirst(==(a),B)
        end
    end
    count
end



#####################################################################
#                          INTERSECTION 3 ELEMENTS
#####################################################################

function find_intersection_indices!(I::Matrix{Int}, A::Int, B::Int, C::Int) 
    if A == B == C
        I[1,1] = 1
        I[2,1] = 1
        I[3,1] = 1
        return 1
    end
    return 0
end

function find_intersection_indices!(I::Matrix{Int}, A::NTuple{N,Int}, B::Int, C::Int) where {N}
    for i ∈ eachindex(A)
        if A[i] == B && A[i] == C
            I[1,1] = i
            I[2,1] = 1
            I[3,1] = 1
            return 1
        end
    end
    return 0
end

function find_intersection_indices!(I::Matrix{Int}, A::Int, B::NTuple{N,Int}, C::Int) where {N}
    for i ∈ eachindex(B)
        if B[i] == A && B[i] == C
            I[1,1] = 1
            I[2,1] = i
            I[3,1] = 1
            return 1
        end
    end
    return 0
end

function find_intersection_indices!(I::Matrix{Int}, A::Int, B::Int, C::NTuple{N,Int}) where {N}
    for i ∈ eachindex(C)
        if C[i] == B && C[i] == A
            I[1,1] = 1
            I[2,1] = 1
            I[3,1] = i
            return 1
        end
    end
    return 0
end


function find_intersection_indices!(I::Matrix{Int}, A::NTuple{N,Int}, B::NTuple{M,Int}, C::Int) where {N,M}
    for i ∈ eachindex(A)
        if A[i] == C && A[i] ∈ B 
            I[1,1] = i
            I[2,1] = findfirst(==(A[i]),B)
            I[3,1] = 1
            return 1
        end
    end
    return 0
end

function find_intersection_indices!(I::Matrix{Int}, A::Int, B::NTuple{M,Int}, C::NTuple{N,Int}) where {N,M}
    for i ∈ eachindex(B)
        if B[i] == A && B[i] ∈ C 
            I[1,1] = 1
            I[2,1] = i
            I[3,1] = findfirst(==(B[i]),C)
            return 1
        end
    end
    return 0
end

function find_intersection_indices!(I::Matrix{Int}, A::NTuple{N,Int}, B::Int, C::NTuple{M,Int}) where {N,M}
    for i ∈ eachindex(A)
        if A[i] == B && A[i] ∈ C 
            I[1,1] = i
            I[2,1] = 1
            I[3,1] = findfirst(==(A[i]),C)
            return 1
        end
    end
    return 0
end

function find_intersection_indices!(I::Matrix{Int}, A::NTuple{N,Int}, B::NTuple{M,Int}, C::NTuple{P,Int}) where {N,M,P}
    count = 0
    for a ∈ A
        if a ∈ B && a ∈ C
            count += 1
            I[1,count] = findfirst(==(a),A)
            I[2,count] = findfirst(==(a),B)
            I[3,count] = findfirst(==(a),C)
        end
    end
    count
end
