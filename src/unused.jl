"""
    get_node_community_weights(node::Int, state::LeidenState, old_comm::Int)

Get weights from node to each neighboring community.
old_comm is the community the node was just removed from.
"""
function get_node_community_weights(node::Int, state::LeidenState, old_comm::Int)
    community_weights = Dict{Int, Float64}()
    
    # Check all neighbors
    for neighbor in neighbors(state.graph, node)
        neighbor_comm = state.partition.membership[neighbor]
        if neighbor_comm == -1  # Skip temporary markers
            continue
        end
        weight = state.weights[node, neighbor]
        
        if haskey(community_weights, neighbor_comm)
            community_weights[neighbor_comm] += weight
        else
            community_weights[neighbor_comm] = weight
        end
    end
    
    # Include old community as an option (even if no edges to it)
    if !haskey(community_weights, old_comm)
        community_weights[old_comm] = 0.0
    end
    
    return community_weights
end