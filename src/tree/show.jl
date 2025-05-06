function Base.Multimedia.display(tree::Root)
    println("A phylogenetic root node (i.e. a tree) with two branches, left and right. The state at the root is \"$(tree.state)\"")
end

function Base.Multimedia.display(node::Node)
    println("An internal node, with three branches, inbounds, left and right. The state at this node is \"$(node.state)\"")
end

function Base.Multimedia.display(tip::Tip)
    println("A tip node with one incoming branch, representing \"$(tip.species_name)\". The state at this tip is \"$(tip.state)\"")
end

function Base.Multimedia.display(branch::Branch)
    println("A branch with a parent node and a descendant node. The states along this branch are $(branch.states), and the state times are $(branch.times), in the order of most recent to most ancient")
end

function Base.Multimedia.display(v::Vector{Union{Root, Node}})
    println("A vector of $(length(v)) internal nodes.")
end

function Base.Multimedia.display(v::Vector{Branch})
    println("A vector of $(length(v)) branches.")
end
