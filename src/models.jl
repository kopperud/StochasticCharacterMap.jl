abstract type Model end

export MkModel
export number_of_states

mutable struct MkModel <: Model
    state_space::Vector{String}
    λ::Float64
    k::Int64
end

function number_of_states(model::Model)
    return(model.k) 
end
