export coordinates
export treeplot

function treeplot(
    tree::Root
)

    fig = CairoMakie.Figure()
    ax = CairoMakie.Axis(fig[1,1], xreversed = true)

    CairoMakie.hidespines!(ax)
    CairoMakie.hidedecorations!(ax)

    ## horizontal lines
    x, y, states = coordinates_horizontal(tree)

    c_states_horizontal = String[]
    points_horizontal = CairoMakie.Point2{Float64}[]
    for index in keys(x)
        for i in 2:length(x[index])
            from = CairoMakie.Point(x[index][i-1], y[index])
            to   = CairoMakie.Point(x[index][i],   y[index])
        
            push!(points_horizontal, from)
            push!(points_horizontal, to)
            push!(c_states_horizontal, states[index][i-1])
        end
    end

    ## vertical lines
    x, y, states = coordinates_vertical(tree)

    c_states_vertical = String[]
    points_vertical = CairoMakie.Point2{Float64}[]
    for index in keys(x)
        #for i in 2:length(x[index])
            from = CairoMakie.Point(x[index], y[index][1])
            to   = CairoMakie.Point(x[index], y[index][2])
        
            push!(points_vertical, from)
            push!(points_vertical, to)
            push!(c_states_vertical, states[index])
        #end
    end
    
    points = vcat(points_horizontal, points_vertical)
    c_states = vcat(c_states_horizontal, c_states_vertical)

    state_space = unique(values(tipstates(tree)))

    tbl = Dict(x => i for (i,x) in enumerate(state_space))
    cs = [tbl[x] for x in c_states]

    cmap = CairoMakie.cgrad([:black, :red, :green], categorical=true)

    CairoMakie.linesegments!(ax, points, color = cs, colormap = cmap, linewidth = 2)

    return(fig)

end



function coordinates_horizontal(tree::Root)
    n_edges = number_of_edges(tree)

    y = zeros(n_edges)
    find_y!(tree.left.outbounds, y)
    find_y!(tree.right.outbounds, y)

    root_height = maximum(node_depths(tree))

    x = Dict{Int64,Vector{Float64}}()
    states = Dict{Int64, Vector{String}}()
    find_x!(tree, x, states, root_height)

    return(x, y, states)
end


function find_y!(node::Tip, y::Vector{Float64})
    index = node.inbounds.index  ## parental branch index 
    y[index] = node.index
end

function find_y!(node::Node, y::Vector{Float64})
    left_branch = node.left
    left_node = left_branch.outbounds

    right_branch = node.right
    right_node = right_branch.outbounds

    find_y!(left_node, y)
    find_y!(right_node, y)

    #index = node.index
    index = node.inbounds.index

    y[index] = (y[left_branch.index] + y[right_branch.index]) / 2
end

function coordinates_vertical(tree::Root)

    nd = node_depths(tree)
    root_height = maximum(nd)

    y = Dict{Int64, Vector{Float64}}()
    x = Dict{Int64, Float64}()
    states = Dict{Int64, String}()

    y_middle = find_vertical!(tree, x, y, states, root_height)

    return(x, y, states)
end

function find_vertical!(node::Tip, x, y, states, t)
    return(node.index)
end

function find_vertical!(node::T, x, y, states, t) where {T <: InternalNode}
    left_branch = node.left
    left_bl = sum(left_branch.times)
    left_node = left_branch.outbounds

    right_branch = node.right
    right_bl = sum(right_branch.times)
    right_node = right_branch.outbounds

    y_left = find_vertical!(left_node, x, y, states, t - left_bl)
    y_right = find_vertical!(right_node, x, y, states, t - right_bl)

    index = node.index

    y[index] = [y_left, y_right]
    x[index] = t
    states[index] = left_branch.states[end]

    y_middle = (y_left + y_right) / 2
    return(y_middle)
end


function find_x!(
    node::T,
    x::Dict{Int64,Vector{Float64}},
    states::Dict{Int64, Vector{String}},
    t::Float64
    ) where {T <: InternalNode}

    for branch in (node.left, node.right)
        bl = sum(branch.times)
        t1 = t - bl
        times = t .- cumsum(reverse(branch.times))
        prepend!(times, t)

        x[branch.index] = times
        states[branch.index] = reverse(branch.states)
        find_x!(branch.outbounds, x, states, t1)
    end
end

function find_x!(::Tip, ::Dict{Int64}, states, t::Float64)
end

export node_depths

## distance from the root to the node
function node_depths(tree::Root)
    n_nodes = number_of_nodes(tree)
    depths = zeros(Float64, n_nodes)

    t = 0.0

    node_depths_po!(tree, depths, t)
    return(depths)
end

function node_depths_po!(node::T, depths::Vector{Float64}, t::Float64) where {T <: InternalNode}
    depths[node.index] = t

    tl = sum(node.left.times)
    node_depths_po!(node.left.outbounds, depths, t + tl)

    tr = sum(node.right.times)
    node_depths_po!(node.right.outbounds, depths, t + tr)
end

function node_depths_po!(node::Tip, depths::Vector{Float64}, t::Float64)
    depths[node.index] = t
end

