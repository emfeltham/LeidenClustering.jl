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

True Leiden refinement (Step B, CPM objective).

For each current community:
  1. Induce the subgraph on that community's nodes.
  2. Initialize a single subcommunity containing all nodes.
  3. Perform local moving within the subgraph using CPM ΔQ:
         ΔQ = [k_i,in(B) - γ*|B|] - [k_i,in(A) - γ*(|A|-1)]
     where A is the node's current subcommunity, B is a candidate target
     (including the option to create a new subcommunity, |B| = 0).
     For weighted graphs, k_i,in(·) sums edge weights; the CPM penalty uses node counts.
  4. Split each resulting subcommunity into connected components (well-connectedness).
  5. Replace the parent community with these subcommunities.
"""
function refine_partition_leiden!(state::LeidenState)
    old_communities = copy(state.partition.community_nodes)

    new_membership = zeros(Int, nv(state.graph))
    new_community_nodes = Dict{Int, Set{Int}}()
    next_id = 1
    γ = state.partition.resolution

    for (_, nodes_in_comm) in old_communities
        # Trivial cases
        if length(nodes_in_comm) <= 1
            new_community_nodes[next_id] = nodes_in_comm
            for v in nodes_in_comm
                new_membership[v] = next_id
            end
            next_id += 1
            continue
        end

        # Build induced subgraph vertex list (global -> local mapping)
        node_list = collect(nodes_in_comm)
        local_index = Dict(node => i for (i, node) in enumerate(node_list))

        # Induced adjacency/weights for the subgraph (dense view over the induced index set).
        # If you keep `state.weights` as SparseMatrixCSC{Float64}, this slice returns a SparseMatrix.
        ws = state.weights[node_list, node_list]

        nloc = length(node_list)

        # --- Subcommunity state (local indices) ---
        # Start with ONE subcommunity that contains all nodes (key = 1)
        sub_membership = fill(1, nloc)   # local node -> subcomm id
        sub_size = Dict(1 => nloc)       # subcomm id -> number of nodes

        # Helper: k_i,in(subcomm) for node v_local under current sub_membership
        k_in_to = function (v_local::Int, target_sc::Int)
            s = 0.0
            # sum weights from v_local to nodes currently in target_sc
            @inbounds for (u, val) in zip(findnz(ws[v_local, :])...)
                # findnz on a row: first tuple element is column index vector, second is values
                # But Graphs/Sparse style: findnz(A[i,:]) returns (J, V)
            end
            # Workaround because zip(findnz(...)) above isn't portable across versions:
            J, V = findnz(ws[v_local, :])
            @inbounds for idx in eachindex(J)
                u = J[idx]
                if sub_membership[u] == target_sc
                    s += V[idx]
                end
            end
            return s
        end

        improved = true
        # Local moving loop inside the community
        while improved
            improved = false

            # Randomized node order (local indices)
            for v_local in shuffle(1:nloc)
                old_sc = sub_membership[v_local]

                # Precompute removal term from old_sc
                k_in_old = k_in_to(v_local, old_sc)
                old_size = sub_size[old_sc]
                # Δ removal part: -(k_i,in(old_sc) - γ*(|old_sc|-1))
                base_removal = -(k_in_old - γ * (old_size - 1))

                # Temporarily remove v_local from old_sc (for size accounting when evaluating targets)
                sub_membership[v_local] = 0
                sub_size[old_sc] = old_size - 1
                if sub_size[old_sc] == 0
                    delete!(sub_size, old_sc)
                end

                # Candidate target subcommunities: neighbors' subcomms + option to create new one
                neighbor_sc = Set{Int}()
                J, _ = findnz(ws[v_local, :])
                @inbounds for u in J
                    sc = sub_membership[u]
                    if sc != 0
                        push!(neighbor_sc, sc)
                    end
                end

                best_sc = old_sc
                best_gain = 0.0

                # 1) Consider moving to each neighboring subcommunity
                for sc in neighbor_sc
                    k_in_new = k_in_to(v_local, sc)
                    size_new = sub_size[sc]  # size BEFORE insertion
                    Δ = base_removal + (k_in_new - γ * size_new)
                    if Δ > best_gain + 1e-12
                        best_gain = Δ
                        best_sc = sc
                    end
                end

                # 2) Consider creating a NEW subcommunity (size_new = 0, k_in_new = 0)
                # Δ_new = base_removal + (0 - γ*0) = base_removal
                if base_removal > best_gain + 1e-12
                    # Create a fresh subcommunity id
                    # Use max key + 1 (avoid relying on contiguity)
                    new_sc = (isempty(sub_size) ? 0 : maximum(keys(sub_size))) + 1
                    best_sc = new_sc
                    best_gain = base_removal
                end

                # Apply the best move (only if strictly improving CPM)
                if best_gain > 1e-12
                    sub_membership[v_local] = best_sc
                    sub_size[best_sc] = get(sub_size, best_sc, 0) + 1
                    if best_sc != old_sc
                        improved = true
                    end
                else
                    # Revert to the old subcommunity
                    sub_membership[v_local] = old_sc
                    sub_size[old_sc] = get(sub_size, old_sc, 0) + 1
                end
            end
        end

        # Group local nodes by subcommunity
        sub_groups = Dict{Int, Set{Int}}()
        for (v_local, sc) in enumerate(sub_membership)
            if sc == 0
                # should not happen; safety
                continue
            end
            if !haskey(sub_groups, sc)
                sub_groups[sc] = Set{Int}()
            end
            push!(sub_groups[sc], v_local)
        end

        # Convert to global ids and enforce connectivity
        for (_, sub_nodes_local) in sub_groups
            # Convert to global vertex ids
            sub_nodes_global = Set( node_list[i] for i in sub_nodes_local )

            # Split into connected components (well-connectedness guarantee)
            for comp in find_connected_components(state.graph, sub_nodes_global)
                new_community_nodes[next_id] = comp
                for v in comp
                    new_membership[v] = next_id
                end
                next_id += 1
            end
        end
    end

    # Commit refined partition
    state.partition.membership = new_membership
    state.partition.community_nodes = new_community_nodes
    rebuild_partition_structures!(state)
    return length(new_community_nodes)
end
