# leiden.jl

"""
    leiden(g::AbstractGraph; resolution::Float64 = 1.0, max_iterations::Int = 100, 
           tol::Float64 = 1e-6, seed::Union{Int, Nothing} = nothing, 
           use_refinement::Bool = true, use_true_leiden::Bool = false)

Full Leiden algorithm with aggregation.
- use_refinement: If true, applies refinement phase (connectivity check at minimum)
- use_true_leiden: If true, uses full Leiden refinement (well-connectedness), otherwise just connectivity
"""
function leiden(g::AbstractGraph; resolution::Float64 = 1.0, max_iterations::Int = 100,
                tol::Float64 = 1e-6, seed::Union{Int, Nothing} = nothing,
                use_refinement::Bool = true, use_true_leiden::Bool = false)
    
    # Keep track of hierarchical structure
    level = 0
    current_graph = g
    current_weights = nothing
    
    # Track membership at each level
    membership_levels = []
    
    prev_modularity = -Inf
    
    while level < max_iterations
        level += 1
        
        # Initialize state for current graph
        state = initialize_state(current_graph, resolution=resolution, weights=current_weights)
        
        # Local moving phase
        local_moving!(state, max_iterations=100, tol=tol, seed=seed)
        
        # Refinement phase (optional)
        if use_refinement
            refine_partition!(state, use_true_leiden=use_true_leiden)
        end
        
        # Check for convergence
        if state.partition.modularity - prev_modularity <= tol
            # Map back to original nodes if we're not at level 1
            if level > 1 && !isempty(membership_levels)
                # Need to map current partition back through all levels
                final_membership = state.partition.membership
                
                # Go through levels in reverse
                for level_membership in reverse(membership_levels)
                    new_membership = zeros(Int, length(level_membership))
                    for (node, comm) in enumerate(level_membership)
                        if comm > 0 && comm <= length(final_membership)
                            new_membership[node] = final_membership[comm]
                        end
                    end
                    final_membership = new_membership
                end
                
                # Renumber communities to be consecutive starting from 1
                unique_comms = unique(final_membership)
                comm_mapping = Dict(old => new for (new, old) in enumerate(sort(unique_comms)))
                final_membership = [comm_mapping[c] for c in final_membership]
                
                # Create final state with original graph
                orig_state = initialize_state(g, resolution=resolution)
                orig_state.partition.membership = final_membership
                
                # Rebuild community structures
                orig_state.partition.community_nodes = Dict{Int, Set{Int}}()
                for (node, comm) in enumerate(final_membership)
                    if !haskey(orig_state.partition.community_nodes, comm)
                        orig_state.partition.community_nodes[comm] = Set{Int}()
                    end
                    push!(orig_state.partition.community_nodes[comm], node)
                end
                
                # CRITICAL: Rebuild weights before computing modularity
                rebuild_partition_structures!(orig_state)
                
                return orig_state
            end
            return state
        end
        
        prev_modularity = state.partition.modularity
        
        # Check if can aggregate further
        if length(state.partition.community_nodes) == nv(current_graph)
            # No change, stop
            break
        end
        
        # Store membership for this level
        push!(membership_levels, copy(state.partition.membership))
        
        # Aggregate to super-graph
        current_graph, current_weights = aggregate_graph(state)
        
        # If super-graph has only one node, we're done
        if nv(current_graph) == 1
            break
        end
    end
    
    # Map memberships back to original nodes if needed
    if level > 1 && !isempty(membership_levels)
        # Get the last state's membership
        if length(state.partition.community_nodes) > 0
            current_membership = collect(1:nv(current_graph))
            for (node, comm) in enumerate(state.partition.membership)
                current_membership[node] = comm
            end
        else
            # If we ended with single community, all nodes in community 1
            current_membership = ones(Int, nv(current_graph))
        end
        
        # Map back through each level
        final_membership = current_membership
        for level_membership in reverse(membership_levels)
            new_membership = zeros(Int, length(level_membership))
            for (node, comm) in enumerate(level_membership)
                if comm > 0 && comm <= length(final_membership)
                    new_membership[node] = final_membership[comm]
                else
                    new_membership[node] = 1  # Default to community 1 if mapping fails
                end
            end
            final_membership = new_membership
        end
        
        # Renumber communities to be consecutive starting from 1
        unique_comms = unique(final_membership)
        comm_mapping = Dict(old => new for (new, old) in enumerate(sort(unique_comms)))
        final_membership = [comm_mapping[c] for c in final_membership]
        
        # Create final state with original graph
        final_state = initialize_state(g, resolution=resolution)
        final_state.partition.membership = final_membership
        
        # Rebuild community structures
        final_state.partition.community_nodes = Dict{Int, Set{Int}}()
        for (node, comm) in enumerate(final_membership)
            if !haskey(final_state.partition.community_nodes, comm)
                final_state.partition.community_nodes[comm] = Set{Int}()
            end
            push!(final_state.partition.community_nodes[comm], node)
        end
        
        # Use the rebuild function to properly set weights and modularity
        rebuild_partition_structures!(final_state)
        
        return final_state
    end
    
    # For level 1, just return the state
    return state
end
