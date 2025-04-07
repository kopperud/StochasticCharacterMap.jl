

export ancestral_character_map

export stochastic_character_map 

function stochastic_character_map(
    tree::Root, 
    model::Model,
    root_prior::Union{Vector{Float64}, Missing} = missing,
    )

    if ismissing(root_prior)
        k = number_of_states(model)
        root_prior = ones(Float64, k) ./ k
    end

    

    #S = ancestral_state_probabilities(tree, model)
    S = ancestral_state_probabilities(model, tree)

    tree2 = deepcopy(tree)

    redraw_nodes_recursive!(tree2, model, S)

    redraw_branches_recursive!(tree2, model)

    return(tree2)
end

