# leiden.jl

"""
    leiden(g::AbstractGraph; resolution::Float64 = 1.0, max_iterations::Int = 100, 
           tol::Float64 = 1e-6, seed::Union{Int, Nothing} = nothing, 
           use_refinement::Bool = true, use_true_leiden::Bool = false)

Full Leiden algorithm with aggregation for community detection.

# Arguments
- `g::AbstractGraph`: The input graph. Can be:
  - `SimpleGraph` for unweighted graphs (all edges have weight 1.0)
  - `SimpleWeightedGraph` for weighted graphs (weights automatically extracted)
- `resolution::Float64 = 1.0`: Resolution parameter (γ) for modularity. Higher values produce more communities.
- `max_iterations::Int = 100`: Maximum number of aggregation levels.
- `tol::Float64 = 1e-6`: Convergence tolerance for modularity improvement.
- `seed::Union{Int, Nothing} = nothing`: Random seed for reproducibility.
- `use_refinement::Bool = true`: If true, applies refinement phase (connectivity check at minimum).
- `use_true_leiden::Bool = false`: If true, uses full Leiden refinement (well-connectedness), 
  otherwise just ensures connectivity.

# Returns
- `LeidenState`: Contains the final partition, communities, and modularity score.

# Notes
- Automatically detects and uses edge weights from `SimpleWeightedGraph` types.
- Edge weights must be non-negative (negative weights will raise an error).
- Zero-weight edges are treated as non-edges.
- The algorithm iteratively aggregates communities until convergence or `max_iterations` is reached.

# Examples
```julia
# Unweighted graph
g = SimpleGraph(100)
# ... add edges ...
state = leiden(g, resolution=1.0)

# Weighted graph  
using SimpleWeightedGraphs
g = SimpleWeightedGraph(100)
add_edge!(g, 1, 2, 5.0)  # edge with weight 5.0
# ... add more weighted edges ...
state = leiden(g, resolution=1.5)  # higher resolution for finer communities
```
"""
function leiden(
    g::AbstractGraph; resolution::Float64 = 1.0,
    max_iterations::Int = 100, tol::Float64 = 1e-6,
    seed::Union{Int,Nothing} = nothing,
    use_refinement::Bool = true, use_true_leiden::Bool = false
)

    # Early exit if no levels requested
    if max_iterations <= 0
        st = initialize_state(g, resolution=resolution)
        return st
    end

    level = 0
    current_graph = g
    current_weights::Union{Nothing,SparseMatrixCSC{Float64}} = nothing
    membership_levels = Vector{Vector{Int}}()
    prev_modularity = -Inf

    last_state = nothing  # keep the most recent state

    while level < max_iterations
        level += 1

        state = initialize_state(current_graph, resolution=resolution, weights=current_weights)
        last_state = state

        # Phase 1
        local_moving!(state, max_iterations=100, tol=tol, seed=seed)
        # Phase 2
        if use_refinement
            refine_partition!(state, use_true_leiden=use_true_leiden)
        end

        # Convergence across levels
        if state.partition.modularity - prev_modularity <= tol
            # If we have previous levels, map back; otherwise just return this state
            if level > 1 && !isempty(membership_levels)
                return _map_back_and_finalize(g, state.partition.membership, membership_levels, resolution)
            else
                return state
            end
        end

        prev_modularity = state.partition.modularity

        # If no change in number of communities, stop
        if length(state.partition.community_nodes) == nv(current_graph)
            break
        end

        # Save membership for this level before aggregation
        push!(membership_levels, copy(state.partition.membership))

        # Phase 3: aggregate
        current_graph, current_weights = aggregate_graph(state)

        if nv(current_graph) == 1
            break
        end
    end

    # Finished due to break or max levels: map back if we built a hierarchy
    if level > 1 && !isempty(membership_levels) && last_state !== nothing
        return _map_back_and_finalize(g, last_state.partition.membership, membership_levels, resolution)
    end

    # Single level case
    return last_state === nothing ? initialize_state(g, resolution=resolution) : last_state
end

function _map_back_and_finalize(
    g::AbstractGraph, top_membership::Vector{Int},
    membership_levels::Vector{Vector{Int}},
    resolution::Float64
)
    # Map memberships down the stack
    final_membership = top_membership
    for lvl in reverse(membership_levels)
        newm = zeros(Int, length(lvl))
        @inbounds for i in eachindex(lvl)
            c = lvl[i]
            newm[i] = (1 <= c <= length(final_membership)) ? final_membership[c] : 1
        end
        final_membership = newm
    end
    # Renumber to 1:N
    uniques = sort(unique(final_membership))
    cmap = Dict(u => i for (i, u) in enumerate(uniques))
    final_membership = [cmap[c] for c in final_membership]

    # Build final state on original graph
    st = initialize_state(g, resolution=resolution)
    st.partition.membership = final_membership
    st.partition.community_nodes = Dict{Int,Set{Int}}()
    for (i, c) in enumerate(final_membership)
        if !haskey(st.partition.community_nodes, c)
            st.partition.community_nodes[c] = Set{Int}()
        end
        push!(st.partition.community_nodes[c], i)
    end
    rebuild_partition_structures!(st)
    return st
end
