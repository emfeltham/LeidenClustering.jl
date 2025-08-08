# local_moving.jl

"""
    local_moving!(state::LeidenState; max_iterations::Int = 100, tol::Float64 = 1e-6, seed::Union{Int, Nothing} = nothing)

Phase 1: Local moving phase of Leiden algorithm.
Move nodes between communities to maximize modularity.
"""
function local_moving!(state::LeidenState; max_iterations::Int = 100, tol::Float64 = 1e-6, seed::Union{Int, Nothing} = nothing)
    if seed !== nothing
        Random.seed!(seed)
    end
    
    n = nv(state.graph)
    improved = true
    iteration = 0
    
    while improved && iteration < max_iterations
        improved = false
        iteration += 1
        
        # Process nodes in random order
        node_order = shuffle(1:n)
        
        for node in node_order
            # Remove node from its community
            old_comm = remove_node_from_community!(node, state)
            
            # Get neighboring communities (pass old_comm to avoid -1 issue)
            comm_weights = get_node_community_weights(node, state, old_comm)
            
            # Find best community
            best_comm = old_comm
            best_gain = 0.0
            
            # Evaluate gain for each neighboring community
            for (comm, _) in comm_weights
                gain = modularity_gain(node, comm, state)
                
                if gain > best_gain + tol
                    best_gain = gain
                    best_comm = comm
                end
            end
            
            # Insert node into best community (could be the same one)
            insert_node_into_community!(node, best_comm, state)
            
            if best_comm != old_comm
                improved = true
            end
        end
        
        # Update modularity
        state.partition.modularity = compute_modularity(state)
    end
    
    return iteration
end

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
