module StochasticCharacterMap

import CairoMakie
import Distributions
import Random

# Write your package code here.


##

## tree data structure
include("tree/types.jl")

##
include("models.jl")

## rest oftree
include("tree/Tree.jl")
include("tree/reindex.jl")
include("tree/number_of_nodes.jl")
include("tree/copy.jl")
include("tree/number_state_changes.jl")
include("tree/treelength.jl")
include("tree/show.jl")
include("tree/find_node.jl")
include("tree/parent_nodes.jl")

## tree moves
include("moves/sample_branch.jl")
include("moves/redraw_recursive.jl")
include("moves/redraw_node.jl")
include("moves/redraw_branch.jl")
include("moves/reassign_copies.jl")

## draw scm
include("moves/stochastic_character_map.jl")

## io
include("readsimmap/readnewick.jl")

## plot
include("plot.jl")

#
include("mk_probs.jl")

end
