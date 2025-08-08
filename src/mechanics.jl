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
Automatically detects and uses weights from SimpleWeightedGraph types.
"""
function initialize_state(g::AbstractGraph; resolution::Float64 = 1.0, wgt::Union{Nothing, SparseMatrixCSC{Float64}} = nothing)
    n = nv(g)
    
    # Create weight matrix - check for SimpleWeightedGraph
    if wgt === nothing
        if isdefined(Main, :SimpleWeightedGraph) && g isa Main.SimpleWeightedGraph
            wgt = weights(g)  # Get weights from the graph
        elseif isdefined(Main, :SimpleWeightedDiGraph) && g isa Main.SimpleWeightedDiGraph
            wgt = weights(g)  # Get weights from the graph
        else
            wgt = sparse(Float64.(adjacency_matrix(g)))  # Default to 1.0 for unweighted
        end
    end
    
    # Validate weights (check for negative values)
    if any(x -> x < 0, nonzeros(wgt))
        throw(ArgumentError("Negative edge weights are not supported in modularity optimization"))
    end
    
    # Calculate node degrees (sum of incident edge weights)
    node_degrees = vec(sum(wgt, dims=2))
    
    # Total weight is sum of all edges (counting each edge once)
    total_weight = sum(wgt) / 2.0
    
    # Initialize partition
    partition = initialize_partition(g, wgt, resolution)
    
    # Initialize community weights
    community_weights_in = Dict{Int, Float64}()
    community_weights_tot = Dict{Int, Float64}()
    
    for i in 1:n
        # Self-loops contribute to internal weight
        community_weights_in[i] = wgt[i, i] / 2.0
        community_weights_tot[i] = node_degrees[i]
    end
    
    state = LeidenState(
        g, wgt, node_degrees, total_weight,
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
    aggregate_graph(state::LeidenState)

Phase 3: Create super-graph where each community becomes a node.
Returns new graph and weights matrix.
If SimpleWeightedGraphs is available, returns a SimpleWeightedGraph.
"""
function aggregate_graph(state::LeidenState)
    comms = sort(collect(keys(state.partition.community_nodes)))
    num_comms = length(comms)
    
    # Calculate weights between communities first
    comm_weights = zeros(num_comms, num_comms)
    
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
            end
            
            # Set weights in matrix (symmetric)
            if weight > 0
                comm_weights[i, j] = weight
                comm_weights[j, i] = weight
            end
        end
    end
    
    # Create the appropriate graph type based on what's available
    if isdefined(Main, :SimpleWeightedGraph)
        # Create a SimpleWeightedGraph
        new_graph = Main.SimpleWeightedGraph(num_comms)
        
        for i in 1:num_comms
            for j in i:num_comms
                if comm_weights[i, j] > 0
                    if i == j
                        # Add self-loop with the weight
                        add_edge!(new_graph, i, j, comm_weights[i, j])
                    else
                        # Add edge with weight (only add once for undirected)
                        add_edge!(new_graph, i, j, comm_weights[i, j])
                    end
                end
            end
        end
        
        # Return the graph and its weight matrix
        new_weights = sparse(comm_weights)
    else
        # Fallback to SimpleGraph if SimpleWeightedGraphs not available
        new_graph = SimpleGraph(num_comms)
        
        for i in 1:num_comms
            for j in i:num_comms
                if i != j && comm_weights[i, j] > 0
                    add_edge!(new_graph, i, j)
                end
            end
        end
        
        new_weights = sparse(comm_weights)
    end
    
    return new_graph, new_weights
end
