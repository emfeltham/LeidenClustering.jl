# mechanics.jl

"""
    initialize_partition(g::AbstractGraph, weights::SparseMatrixCSC{Float64}, resolution::Float64)

Initialize with singleton communities (each node in its own community).
"""
function initialize_partition(g::AbstractGraph, weights::SparseMatrixCSC{Float64}, resolution::Float64)
    n = nv(g)
    membership = collect(1:n)
    community_nodes = Dict{Int, Set{Int}}()
    
    for i in 1:n
        community_nodes[i] = Set([i])
    end
    
    LeidenPartition(membership, community_nodes, -Inf, resolution)
end

"""
    initialize_state(g::AbstractGraph; resolution::Float64 = 1.0, weights::Union{Nothing, SparseMatrixCSC{Float64}} = nothing)

Initialize the Leiden algorithm state.
"""
function initialize_state(g::AbstractGraph; resolution::Float64 = 1.0, weights::Union{Nothing, SparseMatrixCSC{Float64}} = nothing)
    n = nv(g)
    
    # Create weight matrix (default to 1.0 for unweighted graphs)
    if weights === nothing
        weights = sparse(Float64.(adjacency_matrix(g)))
    end
    
    # Calculate node degrees (sum of incident edge weights)
    node_degrees = vec(sum(weights, dims=2))
    
    # Total weight is sum of all edges (counting each edge once)
    total_weight = sum(weights) / 2.0
    
    # Initialize partition
    partition = initialize_partition(g, weights, resolution)
    
    # Initialize community weights
    community_weights_in = Dict{Int, Float64}()
    community_weights_tot = Dict{Int, Float64}()
    
    for i in 1:n
        # Self-loops contribute to internal weight
        community_weights_in[i] = weights[i, i] / 2.0
        community_weights_tot[i] = node_degrees[i]
    end
    
    state = LeidenState(
        g, weights, node_degrees, total_weight,
        partition, community_weights_in, community_weights_tot
    )
    
    # Calculate initial modularity
    state.partition.modularity = compute_modularity(state)
    
    return state
end

"""
    compute_modularity(state::LeidenState)

Compute the modularity of the current partition.
Q = Σ_c [E_c/m - γ*(a_c/(2m))²]
where E_c is internal weight of community c (counted once),
a_c is sum of degrees in c, m is total edge weight, and γ is the resolution parameter.
"""
function compute_modularity(state::LeidenState)
    Q = 0.0
    γ = state.partition.resolution
    m = state.total_weight  # This is already sum(weights)/2
    
    for (comm, E_c) in state.community_weights_in
        a_c = state.community_weights_tot[comm]
        # Correct formula: E_c/m (not E_c/(2m)) because E_c counts internal edges once
        Q += E_c / m
        Q -= γ * (a_c / (2.0 * m))^2
    end
    
    return Q
end

"""
    rebuild_partition_structures!(state::LeidenState)

Rebuild community weight structures based on current partition.
Must be called after updating membership and community_nodes.
"""
function rebuild_partition_structures!(state::LeidenState)
    state.community_weights_in = Dict{Int, Float64}()
    state.community_weights_tot = Dict{Int, Float64}()
    
    for (comm, nodes) in state.partition.community_nodes
        internal_weight = 0.0
        total_weight = 0.0
        
        for node in nodes
            total_weight += state.node_degrees[node]
            
            # Count edges within community
            for neighbor in neighbors(state.graph, node)
                if state.partition.membership[neighbor] == comm
                    internal_weight += state.weights[node, neighbor]
                end
            end
            
            # Add self-loop contribution (diagonal elements)
            internal_weight += state.weights[node, node]
        end
        
        # Divide by 2 because we counted each internal edge twice (once from each endpoint)
        # Self-loops are counted once per node, which is correct
        state.community_weights_in[comm] = internal_weight / 2.0
        state.community_weights_tot[comm] = total_weight
    end
    
    state.partition.modularity = compute_modularity(state)
    return state
end

"""
    remove_node_from_community!(node::Int, state::LeidenState)

Temporarily remove a node from its community (for evaluating moves).
Returns the community it was removed from.
"""
function remove_node_from_community!(node::Int, state::LeidenState)
    comm = state.partition.membership[node]
    
    # Update community weights
    node_degree = state.node_degrees[node]
    state.community_weights_tot[comm] -= node_degree
    
    # Update internal weight - remove edges to other nodes in community
    for neighbor in neighbors(state.graph, node)
        if state.partition.membership[neighbor] == comm && neighbor != node
            state.community_weights_in[comm] -= state.weights[node, neighbor]
        end
    end
    
    # CRITICAL: Remove self-loop contribution
    state.community_weights_in[comm] -= state.weights[node, node] / 2.0
    
    # Remove from community_nodes
    delete!(state.partition.community_nodes[comm], node)
    if isempty(state.partition.community_nodes[comm])
        delete!(state.partition.community_nodes, comm)
        delete!(state.community_weights_in, comm)
        delete!(state.community_weights_tot, comm)
    end
    
    # Mark node as singleton
    state.partition.membership[node] = -1  # Temporary marker
    
    return comm
end

"""
    insert_node_into_community!(node::Int, comm::Int, state::LeidenState)

Insert a node into a community.
"""
function insert_node_into_community!(node::Int, comm::Int, state::LeidenState)
    # Update membership
    state.partition.membership[node] = comm
    
    # Add to community_nodes
    if !haskey(state.partition.community_nodes, comm)
        state.partition.community_nodes[comm] = Set{Int}()
        state.community_weights_in[comm] = 0.0
        state.community_weights_tot[comm] = 0.0
    end
    push!(state.partition.community_nodes[comm], node)
    
    # Update community weights
    node_degree = state.node_degrees[node]
    state.community_weights_tot[comm] += node_degree
    
    # Update internal weight - add edges to other nodes in community
    for neighbor in neighbors(state.graph, node)
        if state.partition.membership[neighbor] == comm
            state.community_weights_in[comm] += state.weights[node, neighbor]
        end
    end
    
    # CRITICAL: Add self-loop contribution
    state.community_weights_in[comm] += state.weights[node, node] / 2.0
end

"""
    modularity_gain(node::Int, comm::Int, state::LeidenState)

Calculate modularity gain from inserting isolated node into community comm.
Node must already be removed from its community.
"""
function modularity_gain(node::Int, comm::Int, state::LeidenState)
    # Weight from node to community
    k_i_in = 0.0
    for neighbor in neighbors(state.graph, node)
        if state.partition.membership[neighbor] == comm
            k_i_in += state.weights[node, neighbor]
        end
    end
    
    # Node degree and community total
    k_i = state.node_degrees[node]
    sigma_tot = get(state.community_weights_tot, comm, 0.0)
    
    # Modularity gain formula
    m = state.total_weight
    resolution = state.partition.resolution
    
    delta_q = k_i_in / m - resolution * (k_i * sigma_tot) / (2.0 * m^2)
    
    return delta_q
end

"""
    refine_partition!(state::LeidenState; use_true_leiden::Bool = false)

Phase 2: Refinement phase.
If use_true_leiden=false: Simple connectivity-based split (Louvain-like)
If use_true_leiden=true: Full Leiden refinement with well-connectedness
"""
function refine_partition!(state::LeidenState; use_true_leiden::Bool = false)
    if !use_true_leiden
        # Simple refinement: just ensure communities are connected
        return refine_partition_simple!(state)
    else
        # True Leiden refinement: ensure well-connectedness
        return refine_partition_leiden!(state)
    end
end

"""
    refine_partition_simple!(state::LeidenState)

Simple refinement - just ensure communities are connected by splitting disconnected ones.
"""
function refine_partition_simple!(state::LeidenState)
    communities = copy(state.partition.community_nodes)
    
    # Clear current partition
    new_membership = zeros(Int, nv(state.graph))
    new_community_nodes = Dict{Int, Set{Int}}()
    new_comm_id = 1
    
    # Check each community for connectivity
    for (comm, nodes) in communities
        if isempty(nodes)
            continue
        end
        
        # If disconnected, split into components
        if !is_connected_community(state.graph, nodes)
            components = find_connected_components(state.graph, nodes)
            for component in components
                new_community_nodes[new_comm_id] = component
                for node in component
                    new_membership[node] = new_comm_id
                end
                new_comm_id += 1
            end
        else
            # Keep as is
            new_community_nodes[new_comm_id] = nodes
            for node in nodes
                new_membership[node] = new_comm_id
            end
            new_comm_id += 1
        end
    end
    
    # Update state
    state.partition.membership = new_membership
    state.partition.community_nodes = new_community_nodes
    
    # Rebuild community weights using the helper function
    rebuild_partition_structures!(state)
    
    return length(new_community_nodes)
end

# Helper: gather weights from node -> neighboring subcommunities (within a parent community)
function _subcomm_weights(node::Int, nodes_in_comm::Set{Int},
                          sub_membership::Dict{Int,Int}, state::LeidenState)
    w = Dict{Int,Float64}()
    for nb in neighbors(state.graph, node)
        if nb in nodes_in_comm
            sc = sub_membership[nb]
            w[sc] = get(w, sc, 0.0) + state.weights[node, nb]
        end
    end
    return w
end

"""
    refine_partition_leiden!(state::LeidenState)

True Leiden refinement on each current community:
- run CPM local moving on the induced subgraph,
- then split subcommunities into connected components,
- replace the parent community by the resulting well-connected subcommunities.
"""
function refine_partition_leiden!(state::LeidenState)
    old_communities = copy(state.partition.community_nodes)

    new_membership = zeros(Int, nv(state.graph))
    new_community_nodes = Dict{Int, Set{Int}}()
    next_id = 1
    γ = state.partition.resolution

    for (_, nodes_in_comm) in old_communities
        # trivial case
        if length(nodes_in_comm) <= 1
            new_community_nodes[next_id] = nodes_in_comm
            for v in nodes_in_comm
                new_membership[v] = next_id
            end
            next_id += 1
            continue
        end

        # --- Initialize singleton subcommunities
        sub_membership = Dict{Int,Int}()              # node -> subcomm id
        sub_sizes      = Dict{Int,Int}()              # subcomm id -> size
        sid = 1
        for v in nodes_in_comm
            sub_membership[v] = sid
            sub_sizes[sid] = 1
            sid += 1
        end

        # --- CPM local moving restricted to nodes_in_comm
        improved = true
        while improved
            improved = false
            for v in shuffle(collect(nodes_in_comm))
                old_sc = sub_membership[v]

                # temporarily remove v from its subcommunity
                sub_sizes[old_sc] -= 1
                if sub_sizes[old_sc] == 0
                    delete!(sub_sizes, old_sc)
                end
                sub_membership[v] = 0  # mark unassigned

                # candidate target subcommunities = neighbor subcomms within this parent community
                neigh_w = _subcomm_weights(v, nodes_in_comm, sub_membership, state)

                best_sc   = old_sc
                best_gain = 0.0

                # Evaluate CPM gain: ΔQ = k_i,in(C) - γ * |C|
                for (sc, w_in) in neigh_w
                    sz = sub_sizes[sc]  # size before insertion
                    Δ = w_in - γ * sz
                    if Δ > best_gain + 1e-12
                        best_gain = Δ
                        best_sc = sc
                    end
                end

                # Option to stay solo (only if better than best_sc)
                # CPM gain to create a new singleton is 0.0 (since |C|=0 and k_i,in=0).
                # So we only move if best_gain > 0.
                if best_gain > 0.0 + 1e-12
                    # insert into best_sc
                    sub_membership[v] = best_sc
                    sub_sizes[best_sc] = get(sub_sizes, best_sc, 0) + 1
                    if old_sc != best_sc
                        improved = true
                    end
                else
                    # revert to original subcommunity
                    if old_sc in keys(sub_sizes)
                        sub_membership[v] = old_sc
                        sub_sizes[old_sc] += 1
                    else
                        # if old_sc got emptied and removed, create a new one for v
                        new_sc = maximum(collect(keys(sub_sizes)); init=0) + 1
                        sub_membership[v] = new_sc
                        sub_sizes[new_sc] = 1
                    end
                end
            end
        end

        # --- Build subcommunities
        sub_groups = Dict{Int, Set{Int}}()
        for (v, sc) in sub_membership
            sub_groups[sc] = get(sub_groups, sc, Set{Int}())
            push!(sub_groups[sc], v)
        end

        # --- Guarantee connectivity of each subcommunity
        for (_, sub_nodes) in sub_groups
            for comp in find_connected_components(state.graph, sub_nodes)
                new_community_nodes[next_id] = comp
                for v in comp
                    new_membership[v] = next_id
                end
                next_id += 1
            end
        end
    end

    # Update state with refined partition and rebuild weights/modularity
    state.partition.membership = new_membership
    state.partition.community_nodes = new_community_nodes
    rebuild_partition_structures!(state)
    return length(new_community_nodes)
end

"""
    aggregate_graph(state::LeidenState)

Phase 3: Create super-graph where each community becomes a node.
Returns new graph, weights matrix, and mapping from new nodes to old communities.
"""
function aggregate_graph(state::LeidenState)
    comms = sort(collect(keys(state.partition.community_nodes)))
    num_comms = length(comms)
        
    # Create new graph
    new_graph = SimpleGraph(num_comms)
    
    # Create weight matrix for new graph
    new_weights = sparse(zeros(num_comms, num_comms))
    
    # Calculate weights between communities
    for i in 1:length(comms)
        for j in i:length(comms)
            comm_i = comms[i]
            comm_j = comms[j]
            
            weight = 0.0
            if i == j
                # Self-loop weight equals twice the internal weight of the community
                # This preserves the internal structure when contracted
                weight = state.community_weights_in[comm_i] * 2.0
            else
                # Inter-community edges: sum all edges between the two communities
                for node_i in state.partition.community_nodes[comm_i]
                    for neighbor in neighbors(state.graph, node_i)
                        if state.partition.membership[neighbor] == comm_j
                            weight += state.weights[node_i, neighbor]
                        end
                    end
                end
                
                # Add edge to new graph if there's a connection
                if weight > 0
                    add_edge!(new_graph, i, j)
                end
            end
            
            # Set weights in matrix (symmetric)
            if weight > 0
                new_weights[i, j] = weight
                new_weights[j, i] = weight
            end
        end
    end
    
    return new_graph, new_weights
end
