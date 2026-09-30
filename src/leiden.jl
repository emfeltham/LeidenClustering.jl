# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Eric Martin Feltham
# Derived from igraph's leiden.c, Copyright (C) 2007-2012 Gabor Csardi
# (GPL-2.0-or-later). See NOTICE.md.

# Leiden algorithm (Traag, Waltman & van Eck 2019), after igraph's leiden.c.
#
# Each pass alternates fast local moving, refinement and aggregation. Intra-cluster weight
# disappears on aggregation, so everything igraph derives from the *original* graph (2m)
# is folded into the resolution once, up front:
#   modularity: node_weights = strength, effective resolution = resolution / 2m
#   CPM:        node_weights = 1,        effective resolution = resolution

# ---------------------------------------------------------------------------
# Phase 1: fast local moving
# ---------------------------------------------------------------------------

# Queue-based local moving. `resolution` is the effective resolution. `membership`
# (ids ≤ n) is updated in place and reindexed to 1:k. Returns (changed, k).
function move_nodes!(
    A::SparseMatrixCSC{Float64,Int},
    membership::Vector{Int},
    node_weights::Vector{Float64},
    resolution::Float64,
    rng::AbstractRNG
)
    n = size(A, 1)
    rows = rowvals(A)
    vals = nonzeros(A)
    changed = false

    cluster_weights = zeros(Float64, n)
    nb_nodes = zeros(Int, n)
    @inbounds for v in 1:n
        c = membership[v]
        cluster_weights[c] += node_weights[v]
        nb_nodes[c] += 1
    end
    empty_clusters = Int[c for c in 1:n if nb_nodes[c] == 0]

    # Circular queue; a node is queued iff it is not stable, so n slots suffice.
    queue = shuffle!(rng, collect(1:n))
    head = 1
    qlen = n
    stable = fill(false, n)
    acc = NeighborAccumulator(n)

    @inbounds while qlen > 0
        v = queue[head]
        head = head == n ? 1 : head + 1
        qlen -= 1

        cur = membership[v]
        nwv = node_weights[v]

        cluster_weights[cur] -= nwv
        nb_nodes[cur] -= 1
        nb_nodes[cur] == 0 && push!(empty_clusters, cur)

        # An empty cluster is always an option
        add!(acc, empty_clusters[end], 0.0)
        for k in nzrange(A, v)
            add!(acc, membership[rows[k]], vals[k])
        end

        best = cur
        max_diff = acc.weight[cur] - nwv * cluster_weights[cur] * resolution
        for c in active(acc)
            diff = acc.weight[c] - nwv * cluster_weights[c] * resolution
            if diff > max_diff
                best = c
                max_diff = diff
            end
        end
        reset!(acc)

        cluster_weights[best] += nwv
        nb_nodes[best] += 1
        best == empty_clusters[end] && pop!(empty_clusters)

        stable[v] = true

        if best != cur
            changed = true
            membership[v] = best
            for k in nzrange(A, v)
                u = rows[k]
                if stable[u] && membership[u] != best
                    stable[u] = false
                    tail = head + qlen
                    tail > n && (tail -= n)
                    queue[tail] = u
                    qlen += 1
                end
            end
        end
    end

    return changed, renumber!(membership)
end

# ---------------------------------------------------------------------------
# Phase 2: refinement
# ---------------------------------------------------------------------------

# Refine every cluster of `membership` into well-connected sub-clusters, starting from
# singletons and randomly merging (probability ∝ exp(gain/β)). Sub-clusters never cross
# cluster boundaries. Returns (refined, k) with ids 1:k.
function refine_partition(
    A::SparseMatrixCSC{Float64,Int},
    membership::Vector{Int},
    node_weights::Vector{Float64},
    resolution::Float64,
    beta::Float64,
    rng::AbstractRNG
)
    n = size(A, 1)
    rows = rowvals(A)
    vals = nonzeros(A)

    nb_clusters = maximum(membership)
    start, members = group_members(membership, nb_clusters)

    # Scratch shared across clusters; entries are indexed by node id, which doubles as
    # the id of the initial singleton refined cluster.
    refined = zeros(Int, n)
    out = zeros(Int, n)
    remap = zeros(Int, n)
    cluster_weights = zeros(Float64, n)
    ext_weight = zeros(Float64, n)
    non_singleton = fill(false, n)
    cum_trans = zeros(Float64, n)
    order = Vector{Int}(undef, n)
    acc = NeighborAccumulator(n)
    nb_refined = 0

    for comm in 1:nb_clusters
        lo, hi = start[comm], start[comm + 1] - 1
        size_c = hi - lo + 1
        if size_c == 1
            nb_refined += 1
            out[members[lo]] = nb_refined
            continue
        end

        subset = view(members, lo:hi)
        total_nw = 0.0
        @inbounds for v in subset
            refined[v] = v
            cluster_weights[v] = node_weights[v]
            non_singleton[v] = false
            total_nw += node_weights[v]
            s = 0.0
            for k in nzrange(A, v)
                membership[rows[k]] == comm && (s += vals[k])
            end
            ext_weight[v] = s
        end

        shuffled = view(order, 1:size_c)
        copyto!(shuffled, subset)
        shuffle!(rng, shuffled)
        @inbounds for v in shuffled
            cur = refined[v]
            if !non_singleton[cur] &&
               ext_weight[cur] >= cluster_weights[cur] * (total_nw - cluster_weights[cur]) * resolution

                nwv = node_weights[v]
                cluster_weights[cur] = 0.0

                add!(acc, cur, 0.0)
                for k in nzrange(A, v)
                    u = rows[k]
                    membership[u] == comm && add!(acc, refined[u], vals[k])
                end

                best = cur
                max_diff = 0.0
                total_cum = 0.0
                ids = active(acc)
                for (j, c) in enumerate(ids)
                    if ext_weight[c] >= cluster_weights[c] * (total_nw - cluster_weights[c]) * resolution
                        diff = acc.weight[c] - nwv * cluster_weights[c] * resolution
                        if diff > max_diff
                            best = c
                            max_diff = diff
                        end
                        diff >= 0 && beta > 0 && (total_cum += exp(diff / beta))
                    end
                    cum_trans[j] = total_cum
                end

                chosen = best
                if beta > 0 && total_cum < Inf
                    r = rand(rng) * total_cum
                    idx = searchsortedfirst(view(cum_trans, 1:length(ids)), r)
                    chosen = ids[min(idx, length(ids))]
                end
                reset!(acc)

                cluster_weights[chosen] += nwv
                for k in nzrange(A, v)
                    u = rows[k]
                    membership[u] == comm || continue
                    if refined[u] == chosen
                        ext_weight[chosen] -= vals[k]
                    else
                        ext_weight[chosen] += vals[k]
                    end
                end

                if chosen != cur
                    refined[v] = chosen
                    non_singleton[chosen] = true
                end
            end
        end

        # Renumber this cluster's refined ids into the global 1:nb_refined range
        @inbounds for v in subset
            c = refined[v]
            if remap[c] == 0
                nb_refined += 1
                remap[c] = nb_refined
            end
            out[v] = remap[c]
        end
        @inbounds for v in subset
            remap[refined[v]] = 0
        end
    end

    return out, nb_refined
end

# ---------------------------------------------------------------------------
# Phase 3: aggregation
# ---------------------------------------------------------------------------

# Contract each refined cluster to a node, dropping intra-cluster weight. The aggregate
# node inherits the (non-refined) cluster of its members, which seeds the next round of
# local moving.
function aggregate(
    A::SparseMatrixCSC{Float64,Int},
    refined::Vector{Int},
    nb_refined::Int,
    node_weights::Vector{Float64},
    membership::Vector{Int}
)
    agg_A, _ = contract(A, refined, nb_refined)
    agg_nw = zeros(Float64, nb_refined)
    agg_memb = zeros(Int, nb_refined)
    @inbounds for v in eachindex(refined)
        c = refined[v]
        agg_nw[c] += node_weights[v]
        agg_memb[c] = membership[v]
    end
    return agg_A, agg_nw, agg_memb
end

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

# One full Leiden run to convergence, starting from `membership0`.
# Returns (membership, changed).
function leiden_iteration(
    A0::SparseMatrixCSC{Float64,Int},
    node_weights0::Vector{Float64},
    membership0::Vector{Int},
    resolution::Float64,
    beta::Float64,
    rng::AbstractRNG
)
    n = size(A0, 1)
    A = A0
    nw = node_weights0
    memb = copy(membership0)
    renumber!(memb)
    result = copy(memb)
    aggregate_node = collect(1:n)   # original node -> node at the current level
    level = 0
    changed_any = false

    while true
        changed, nb_clusters = move_nodes!(A, memb, nw, resolution, rng)
        changed_any |= changed

        # Level 0 moves the original nodes, so its result stands even if we stop here
        # (igraph moves them in place)
        level == 0 && copyto!(result, memb)

        # Stop when every cluster is a single node: nothing left to aggregate
        nb_clusters < size(A, 1) || break

        if level > 0
            @inbounds for i in 1:n
                result[i] = memb[aggregate_node[i]]
            end
        end

        refined, nb_refined = refine_partition(A, memb, nw, resolution, beta, rng)

        # If refinement merged nothing, aggregate on the actual clustering
        if nb_refined >= size(A, 1)
            refined = copy(memb)
            nb_refined = nb_clusters
        end

        @inbounds for i in 1:n
            aggregate_node[i] = refined[aggregate_node[i]]
        end

        A, nw, memb = aggregate(A, refined, nb_refined, nw, memb)
        level += 1
    end

    return result, changed_any
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    leiden_clustering(graph; resolution=1.0, beta=0.01, n_iterations=2,
                      objective=:modularity, weights=Graphs.weights(graph),
                      node_weights=nothing, initial_membership=nothing,
                      rng=Random.default_rng(), seed=nothing)

Leiden community detection (Traag, Waltman & van Eck 2019), after igraph's
`igraph_community_leiden`. Communities are guaranteed to be connected.

# Arguments
- `graph::AbstractGraph`: undirected graph.
- `resolution`: γ; higher gives more, smaller communities.
- `beta`: randomness of the refinement step (`0` is greedy).
- `n_iterations`: maximum number of full passes, each started from the previous result;
  stops early once a pass changes nothing (R/igraph's default is 2). A negative value runs
  until a pass changes nothing; `0` returns `initial_membership` unchanged.
- `objective`: `:modularity`, or `:cpm` (Constant Potts Model).
- `weights`: edge weights, `weights[i, j]` for edge `i–j` (non-negative). The default uses the
  weights of a `SimpleWeightedGraph` (or any graph type that defines `Graphs.weights`), and 1
  for every edge otherwise.
- `node_weights`: vertex weights `n_i` for `:cpm` (default 1 each). Not allowed with
  `:modularity`, which always uses vertex strengths.
- `initial_membership`: starting partition, one positive community id per vertex (default:
  every vertex alone).
- `rng`: random number generator. `seed`, if given, creates a private `Xoshiro(seed)`
  instead; the global RNG is never reseeded.

# Returns
A [`Partition`](@ref); it destructures as `membership, quality = leiden_clustering(g)`.

# Example
```julia
using Graphs, LeidenClustering
result = leiden_clustering(smallgraph(:karate); seed=42)
result.membership, result.quality
```

# References
Traag, V.A., Waltman, L. & van Eck, N.J. (2019). From Louvain to Leiden: guaranteeing
well-connected communities. Scientific Reports 9, 5233.
"""
function leiden_clustering(
    graph::AbstractGraph;
    resolution::Real=1.0,
    beta::Real=0.01,
    n_iterations::Integer=2,
    objective::Symbol=:modularity,
    weights=Graphs.weights(graph),
    node_weights::Union{AbstractVector{<:Real},Nothing}=nothing,
    initial_membership::Union{AbstractVector{<:Integer},Nothing}=nothing,
    rng::AbstractRNG=Random.default_rng(),
    seed::Union{Integer,Nothing}=nothing
)
    is_directed(graph) && throw(ArgumentError("Leiden algorithm requires undirected graph"))
    resolution ≥ 0 || throw(ArgumentError("Resolution must be non-negative"))
    0 ≤ beta ≤ 1 || throw(ArgumentError("Beta must be in [0, 1]"))
    objective ∈ (:modularity, :cpm) || throw(ArgumentError("objective must be :modularity or :cpm"))
    modularity_objective = objective === :modularity
    n = nv(graph)
    if node_weights !== nothing
        modularity_objective &&
            throw(ArgumentError("node_weights apply to objective=:cpm only; modularity uses strengths"))
        length(node_weights) == n ||
            throw(ArgumentError("node_weights must have one entry per vertex"))
        all(w -> isfinite(w) && w ≥ 0, node_weights) ||
            throw(ArgumentError("node_weights must be non-negative and finite"))
    end
    if initial_membership !== nothing
        length(initial_membership) == n ||
            throw(ArgumentError("initial_membership must have one entry per vertex"))
        all(≥(1), initial_membership) ||
            throw(ArgumentError("initial_membership ids must be positive"))
    end
    rng = resolve_rng(rng, seed)
    γ, β = Float64(resolution), Float64(beta)

    n == 0 && return Partition(Int[], 0.0, Float64[])

    adj = Adjacency(graph, weights)
    strength = strengths(adj)
    nw = modularity_objective ? strength :
         node_weights === nothing ? ones(Float64, n) : Vector{Float64}(node_weights)
    twom = sum(strength)
    γ_eff = modularity_objective ? (twom > 0 ? γ / twom : 0.0) : γ

    membership = initial_membership === nothing ? collect(1:n) : Vector{Int}(initial_membership)
    renumber!(membership)
    iteration = 0
    while n_iterations < 0 || iteration < n_iterations
        membership, changed = leiden_iteration(adj.A, nw, membership, γ_eff, β, rng)
        iteration += 1
        changed || break
    end

    quality = modularity_objective ?
        modularity(adj, membership; resolution=γ) :
        cpm_quality(adj, membership, nw, γ)
    return Partition(membership, quality, [quality])
end
