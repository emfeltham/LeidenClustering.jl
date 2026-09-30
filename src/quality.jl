# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Eric Martin Feltham
# Derived in part from igraph's modularity.c, Copyright (C) 2007-2020 The igraph development
# team (GPL-2.0-or-later). See NOTICE.md.

# Partition quality functions.

"""
    modularity(adj::Adjacency, membership; resolution=1.0)

Newman modularity `Q = 1/2m Σ_ij (A_ij − γ k_i k_j / 2m) δ(c_i, c_j)`. Community ids need
not be consecutive. Returns `0.0` for a graph without edges.
"""
function modularity(
    adj::Adjacency,
    membership::AbstractVector{<:Integer};
    resolution::Real=1.0
)
    A = adj.A
    n = size(A, 1)
    length(membership) == n ||
        throw(ArgumentError("Membership vector length must equal number of vertices"))
    resolution ≥ 0 || throw(ArgumentError("Resolution parameter must be non-negative"))
    n == 0 && return 0.0
    minimum(membership) ≥ 1 || throw(ArgumentError("Community ids must be positive"))

    rows = rowvals(A)
    vals = nonzeros(A)
    total = zeros(Float64, maximum(membership))   # strength per community
    inside = 0.0                                  # weight inside communities, 2× per edge
    @inbounds for j in 1:n
        cj = membership[j]
        s = 2 * adj.loops[j]
        inside += s
        for p in nzrange(A, j)
            s += vals[p]
            membership[rows[p]] == cj && (inside += vals[p])
        end
        total[cj] += s
    end
    twom = sum(total)
    twom == 0 && return 0.0
    return (inside - resolution * sum(abs2, total) / twom) / twom
end

"""
    modularity(graph, membership; resolution=1.0, weights=Graphs.weights(graph))

Modularity of the partition `membership` (one community id per vertex) of an undirected
graph. `weights[i, j]` gives the weight of edge `i–j`; the default uses the edge weights of
a `SimpleWeightedGraph` (or any graph type that defines `Graphs.weights`), and 1 otherwise.
Weights must be non-negative. Not exported, because
`Graphs.modularity` has the same name: call it as `LeidenClustering.modularity`.

# Example
```julia
using Graphs, LeidenClustering
g = smallgraph(:karate)
result = leiden_clustering(g; seed=1)
LeidenClustering.modularity(g, result.membership)
```

Based on Newman (2006), *Modularity and community structure in networks*, PNAS 103, 8577.
"""
function modularity(
    graph::AbstractGraph,
    membership::AbstractVector{<:Integer};
    resolution::Real=1.0,
    weights=Graphs.weights(graph)
)
    is_directed(graph) && throw(ArgumentError("modularity is only defined here for undirected graphs"))
    return modularity(Adjacency(graph, weights), membership; resolution)
end

"""
    cpm_quality(adj, membership, node_weights, resolution)

Constant Potts Model quality `1/2m Σ_ij (A_ij − γ n_i n_j) δ(c_i, c_j)`, as igraph computes
it: a self-loop of weight `w` counts `2w`, like in the strengths.
"""
function cpm_quality(
    adj::Adjacency,
    membership::Vector{Int},
    node_weights::AbstractVector{Float64},
    resolution::Float64
)
    A = adj.A
    rows = rowvals(A)
    vals = nonzeros(A)
    total = 2 * sum(adj.loops)
    quality = total
    @inbounds for v in 1:size(A, 1), p in nzrange(A, v)
        total += vals[p]
        membership[rows[p]] == membership[v] && (quality += vals[p])
    end
    cw = zeros(Float64, maximum(membership))
    @inbounds for i in eachindex(membership)
        cw[membership[i]] += node_weights[i]
    end
    quality -= resolution * sum(abs2, cw)
    return total == 0 ? 0.0 : quality / total
end
