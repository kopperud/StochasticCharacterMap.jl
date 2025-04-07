export transition_probability

sample(rng::Random.AbstractRNG, a::Random.AbstractArray) = a[rand(rng, 1:length(a))]
sample(a::AbstractArray) = sample(Random.default_rng(), a)

function logpdf(tree::Root, model::Model)
   
    λ = model.λ
    k = model.k

    n_branches = number_of_edges(tree)
    D = zeros(eltype(λ), n_branches, 2, k)

    x, log_nf = postorder!(tree, model, D)

    ## what are appropriate root freqs?
    ## the ratio of the observed character states?
    root_freqs = [1/k for _ in 1:k]

    logl = log(sum(root_freqs .* x)) + log_nf
    return(logl)
end


function Base.Multimedia.display(model::MkModel)
    λ = model.λ

    println("A Markov k-state tree distribution with rate λ = $λ.")
end

#######################################
##
##   transition probabilities, i.e. the P(Δt) matrix 
##
########################################
function transition_probability(
    model::Model, 
    t::Float64
    )
    k = model.k
    λ = model.λ

    pii = (1/k) + ((k-1)/k)*exp(-k * λ * t)
    pij = (1/k) - (1/k)*exp(-k * λ * t)

    Q = zeros(eltype(λ), k, k)
    Q[:,:] .= pij
    for i in 1:k
        Q[i,i] = pii
    end
    return(Q)
end


#######################################
##
##               postorder
##
########################################
function postorder!(
    node::N, 
    model::Model,
    D::Array{T, 3}
    ) where {N <: InternalNode, T <: Real}


    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(model, left_bl)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(model, right_bl)

    x_left, log_nf_left = postorder!(left_node, model, D)
    x_right, log_nf_right = postorder!(right_node, model, D)

    D[left_branch.index,1,:] = x_left
    D[right_branch.index,1,:] = x_right
    
    state_left = P_left * x_left
    state_right = P_right * x_right

    D[left_branch.index,2,:] = state_left
    D[right_branch.index,2,:] = state_right
    
    state = state_left .* state_right
    nf_node = sum(state)
    state = state ./ nf_node
    log_nf = log_nf_left + log_nf_right + log(nf_node)

    return(state, log_nf)
end

function postorder!(
        node::Tip, 
        model::Model,
        D::Array{T, 3}
    ) where {T <: Real}

    state = zeros(model.k)
    idx = findfirst(isequal(node.state), model.state_space)
    state[idx] = 1
    log_nf = 0.0

    return(state, log_nf)
end

#######################################
##
##               preorder
##
########################################
function preorder!(
        node::T, 
        model::Model, 
        D::Array{Float64, 3}, 
        F::Array{Float64, 3}, 
        S_parent::Vector{Float64}
    ) where {T <: InternalNode}
    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(model, left_bl)
    F_left = S_parent ./ D[left_branch.index,2,:]
    F[left_branch.index,2,:] = F_left
    F[left_branch.index,1,:] = P_left * F_left
    S_left = F[left_branch.index,1,:] .* D[left_branch.index,1,:]
    S_left = S_left ./ sum(S_left)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(model, right_bl)
    F_right = S_parent ./ D[right_branch.index,2,:]
    F[right_branch.index,2,:] = F_right
    F[right_branch.index,1,:] = P_right * F_right
    S_right = F[right_branch.index,1,:] .* D[right_branch.index,1,:]
    S_right = S_right ./ sum(S_right)

    preorder!(left_node, model, D, F, S_left)
    preorder!(right_node, model, D, F, S_right)
end

function preorder!(
    node::Tip, 
    model::Model, 
    D::Array{Float64, 3}, 
    F::Array{Float64, 3}, 
    S::Vector{Float64}
    )
    ## do nothing
    nothing
end


#######################################
##
##   ancestral state probabilities
##
########################################
export ancestral_state_probabilities

## these are for edge indices
function ancestral_state_probabilities(
    model::Model,
    root::Root,
)
    λ = model.λ
    k = model.k
    n_branches = number_of_edges(root)
    F = zeros(eltype(λ), n_branches, 2, k) ## 2 because beginning and end of branch
    D = zeros(eltype(λ), n_branches, 2, k)

    ## Postorder    
    x, log_nf = postorder!(root, model, D) ## this fills in D

    ## some shenanigans at the root
    left_branch_idx = root.left.index
    right_branch_idx = root.right.index

    F_root = ones(k)
    D_root = D[left_branch_idx,2,:] .* D[right_branch_idx,2,:]
    S_root = F_root .* D_root
    S_root = S_root ./ sum(S_root)

    ## Preorder
    preorder!(root, model, D, F, S_root)
   
    ## ancestral state probs
    S_branches = D[:,1,:] .* F[:,1,:]

    ####################################
    ##
    ##   convert to node indices
    ##
    ##########################
    n_nodes = number_of_nodes(root)
    S_nodes = zeros(n_nodes, k)

    S_nodes[root.index,:] .= S_root
    asp_po(root, S_branches, S_nodes)

    return(S_nodes)
end

function asp_po(
        node::T, 
        S_branches::Matrix{Float64}, 
        S_nodes::Matrix{Float64}
    ) where {T <: InternalNode}

    left_branch_idx = node.left.index
    left_node_idx = node.left.outbounds.index
    S_nodes[left_node_idx,:] = S_branches[left_branch_idx,:]

    right_branch_idx = node.right.index
    right_node_idx = node.right.outbounds.index
    S_nodes[right_node_idx,:] = S_branches[right_branch_idx,:]

    asp_po(node.left.outbounds, S_branches, S_nodes)
    asp_po(node.right.outbounds, S_branches, S_nodes)

end

function asp_po(
        node::Tip, 
        S_branches::Matrix{Float64}, 
        S_nodes::Matrix{Float64}
    )
    node_idx = node.index
    S_nodes[node_idx,:] = S_branches[node.inbounds.index,:]
end
