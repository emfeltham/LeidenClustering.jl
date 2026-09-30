# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Eric Martin Feltham
# Derived from igraph's louvain.c, Copyright (C) 2007-2020 The igraph development team
# (GPL-2.0-or-later). See NOTICE.md.

# Louvain algorithm (Blondel et al. 2008), after igraph's louvain.c.
#
# Each level runs local moving on the current graph, then contracts communities into
# vertices. Intra-community weight is kept as self-loops, so the modularity bookkeeping
# stays exact across levels.
#
# Deviation from igraph: local moving uses a work queue (revisit a vertex only after a
# neighbour moves) instead of repeated full sweeps. Same strictly-positive-gain criterion,
# far fewer vertex visits on large graphs.

# Modularity from per-community totals (communities with size 0 are skipped).
function modularity_from_totals(csize, w_in, w_all, twom, resolution)
    q = 0.0
    @inbounds for c in eachindex(csize)
        csize[c] > 0 || continue
        q += (w_in[c] - resolution * w_all[c] * w_all[c] / twom) / twom
    end
    return q
end

# One Louvain level. Returns (contracted graph, membership of this level's vertices in 1:k,
# modularity of that membership).
function louvain_level(adj::Adjacency, resolution::Float64, rng::AbstractRNG)
    A = adj.A
    n = size(A, 1)
    rows = rowvals(A)
    vals = nonzeros(A)

    strength = strengths(adj)
    twom = sum(strength)
    membership = collect(1:n)
    csize = ones(Int, n)
    w_all = copy(strength)          # total strength of each community
    w_in = 2 .* adj.loops           # weight inside each community, 2× per edge

    # Work queue instead of igraph's repeated full sweeps: a vertex is revisited only after
    # a neighbour moved. Each move strictly raises modularity, so this terminates at a local
    # optimum. The queue is circular; a vertex is queued at most once, so n slots suffice.
    queue = shuffle!(rng, collect(1:n))
    queued = fill(true, n)
    head = 1
    qlen = n
    acc = NeighborAccumulator(n)

    @inbounds while qlen > 0
        v = queue[head]
        head = head == n ? 1 : head + 1
        qlen -= 1
        queued[v] = false

        old = membership[v]
        s_all = strength[v]
        s_loop = 2 * adj.loops[v]
        for p in nzrange(A, v)
            add!(acc, membership[rows[p]], vals[p])
        end
        s_in = acc.weight[old]

        # Take v out of its community
        csize[old] -= 1
        w_all[old] -= s_all
        w_in[old] -= 2 * s_in + s_loop

        # Move only for a strictly positive gain
        best = old
        max_gain = 0.0
        best_weight = s_in
        for c in active(acc)
            w = acc.weight[c]
            gain = w - resolution * w_all[c] * s_all / twom
            if gain > max_gain
                max_gain = gain
                best = c
                best_weight = w
            end
        end
        reset!(acc)

        membership[v] = best
        csize[best] += 1
        w_all[best] += s_all
        w_in[best] += 2 * best_weight + s_loop

        if best != old
            for p in nzrange(A, v)
                u = rows[p]
                if !queued[u] && membership[u] != best
                    queued[u] = true
                    tail = head + qlen
                    tail > n && (tail -= n)
                    queue[tail] = u
                    qlen += 1
                end
            end
        end
    end
    q = modularity_from_totals(csize, w_in, w_all, twom, resolution)

    k = renumber!(membership)
    B, internal = contract(A, membership, k)
    loops = internal
    @inbounds for v in 1:n
        loops[membership[v]] += adj.loops[v]
    end
    return Adjacency(B, loops), membership, q
end

"""
    louvain_clustering(graph; resolution=1.0, max_iterations=100,
                       weights=Graphs.weights(graph), rng=Random.default_rng(), seed=nothing)

Multilevel Louvain community detection maximising modularity, after igraph's
`igraph_community_multilevel`.

Each level moves vertices greedily between communities, then contracts communities into
single vertices, until no level improves modularity. Communities are not guaranteed to be
connected; see [`leiden_clustering`](@ref) if that matters.

# Arguments
- `graph::AbstractGraph`: undirected graph.
- `resolution`: γ; higher gives more, smaller communities.
- `max_iterations`: maximum number of aggregation levels.
- `weights`: edge weights, `weights[i, j]` for edge `i–j` (non-negative). The default uses the
  weights of a `SimpleWeightedGraph` (or any graph type that defines `Graphs.weights`), and 1
  for every edge otherwise.
- `rng`: random number generator. `seed`, if given, creates a private `Xoshiro(seed)`
  instead; the global RNG is never reseeded.

# Returns
A [`Partition`](@ref); `qualities` holds the modularity after each level (the last entry
belongs to `membership`). It destructures as `membership, quality = louvain_clustering(g)`.

# Example
```julia
using Graphs, LeidenClustering
result = louvain_clustering(smallgraph(:karate); seed=1)
```

# References
Blondel, Guillaume, Lambiotte & Lefebvre (2008), *Fast unfolding of communities in large
networks*, J. Stat. Mech. P10008.
"""
function louvain_clustering(
    graph::AbstractGraph;
    resolution::Real=1.0,
    max_iterations::Integer=100,
    weights=Graphs.weights(graph),
    rng::AbstractRNG=Random.default_rng(),
    seed::Union{Integer,Nothing}=nothing
)
    is_directed(graph) && throw(ArgumentError("Louvain algorithm works for undirected graphs only"))
    resolution ≥ 0 || throw(ArgumentError("Resolution parameter must be non-negative"))
    max_iterations ≥ 1 || throw(ArgumentError("Max iterations must be positive"))
    rng = resolve_rng(rng, seed)
    γ = Float64(resolution)

    adj0 = Adjacency(graph, weights)
    membership = collect(1:nv(adj0))
    qualities = Float64[]

    current = adj0
    q = -1.0
    for _ in 1:max_iterations
        n_step = nv(current)
        prev_q = q
        next, step_membership, q = louvain_level(current, γ, rng)
        (nv(next) == n_step || q < prev_q) && break

        @inbounds for i in eachindex(membership)
            membership[i] = step_membership[membership[i]]
        end
        push!(qualities, q)
        current = next
    end

    quality = isempty(qualities) ? modularity(adj0, membership; resolution=γ) : qualities[end]
    return Partition(membership, quality, qualities)
end
