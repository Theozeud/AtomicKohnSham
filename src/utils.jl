struct NothingFunction <: Function end
(::NothingFunction)(args...;kwargs...) = nothing

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

function intersection_with_indices(A, B)
    # Find intersection of A and B
    intersect_elements = intersect(A, B)
    # Create a vector of tuples (element, index_in_A, index_in_B)
    return  [(findfirst(x -> x == el, A), findfirst(x -> x == el, B)) for el in intersect_elements]
end

function intersection_with_indices(A, B, C)
    # Find intersection of A and B
    intersect_elements = intersect(A, B, C)
    # Create a vector of tuples (element, index_in_A, index_in_B)
    return  [(findfirst(x -> x == el, A), findfirst(x -> x == el, B), findfirst(x -> x == el, C)) for el in intersect_elements]
end


function findfirsttwo(predicate, P)
    result = Int[]
    for (i, p) in pairs(P)
        predicate(p) && push!(result, i)
        length(result) == 2 && break
    end
    return result
end