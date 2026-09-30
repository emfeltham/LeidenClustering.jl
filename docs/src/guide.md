# Guide

```@setup guide
using Graphs, SimpleWeightedGraphs, LeidenClustering, Random
```

## Choosing an algorithm

| | [`leiden_clustering`](@ref) | [`louvain_clustering`](@ref) |
|---|---|---|
| Communities connected | guaranteed | not guaranteed |
| Objectives | modularity, CPM | modularity |
| Speed | fast | faster (roughly 2× on large graphs) |
| Result | one final partition | one final partition, plus the modularity of every level |

Use Leiden by default. Louvain can leave a community internally disconnected, because a
vertex moved out of a community can strand others; Leiden's refinement step prevents this.

## Results

Both functions return a [`Partition`](@ref). Community ids are `1:k`, numbered in order of
first appearance.

```@example guide
g = smallgraph(:karate)
result = leiden_clustering(g; seed=1)
result.membership
```

```@example guide
ncommunities(result), result.quality
```

A `Partition` also destructures like a tuple, `membership, quality = result`. The field
`qualities` holds the modularity after each level for Louvain (a one-element vector for
Leiden):

```@example guide
membership, quality = louvain_clustering(g; seed=1)
louvain_clustering(g; seed=1).qualities
```

The community of vertex `v` is `result.membership[v]`. To list the vertices of each community:

```@example guide
communities = [findall(==(c), result.membership) for c in 1:ncommunities(result)]
```

To score any partition, including one from another method, use
`LeidenClustering.modularity(g, membership)`. It is not exported because `Graphs.modularity`
has the same name.

```@example guide
LeidenClustering.modularity(g, result.membership; resolution=1.0)
```

## Weighted graphs

By default edge weights come from `Graphs.weights(g)`: the weights of a
`SimpleWeightedGraph`, or of any other graph type that defines `Graphs.weights`. Graphs
without weights count each edge as 1. Weights must be non-negative and finite.

```@example guide
w = SimpleWeightedGraph(6)
for (a, b, x) in [(1, 2, 5.0), (2, 3, 5.0), (1, 3, 5.0),
                  (4, 5, 5.0), (5, 6, 5.0), (4, 6, 5.0),
                  (3, 4, 0.5)]
    add_edge!(w, a, b, x)
end
leiden_clustering(w; seed=1).membership
```

The two triangles are heavy and the bridge between them is light, so they form two communities.

To use other weights, pass any matrix with `weights[i, j]` the weight of edge `i–j` (entries
for non-edges are ignored). The same keyword works for `louvain_clustering` and
`LeidenClustering.modularity`:

```@example guide
W = zeros(6, 6)
for e in edges(w)
    W[src(e), dst(e)] = W[dst(e), src(e)] = 1.0
end
leiden_clustering(w; weights=W, seed=1).membership   # the unweighted graph
```

## Resolution

`resolution` (γ) trades community size against number. At γ = 0 nothing penalises large
communities; larger values give more, smaller ones.

```@example guide
g = smallgraph(:karate)
[(γ = γ, k = ncommunities(leiden_clustering(g; resolution=γ, seed=1)))
 for γ in (0.25, 0.5, 1.0, 2.0, 4.0)]
```

## Objectives

`objective=:modularity` (the default) compares each community's internal weight with what a
random graph with the same degrees would give. `objective=:cpm` uses the Constant Potts
Model, which compares it with a fixed density instead: `resolution` is then the minimum
internal density of a community, and it avoids modularity's resolution limit on large graphs.
CPM gives unit weight to every vertex, so useful values of `resolution` are much smaller
than for modularity.

```@example guide
[(γ = γ, k = ncommunities(leiden_clustering(g; objective=:cpm, resolution=γ, seed=1)))
 for γ in (0.05, 0.1, 0.2, 0.4)]
```

`node_weights` gives CPM vertex weights ``n_i`` in place of 1, so the penalty for a community
is ``\gamma (\sum_{i \in c} n_i)^2``. It is not allowed with modularity, whose vertex weights
are always the weighted degrees.

## Reproducibility

Both algorithms are randomised. Two ways to control the randomness:

- `seed=42` runs with a private `Xoshiro(42)`. It is fully reproducible for a given Julia
  version and never touches the global random number generator.
- `rng=my_rng` uses any `AbstractRNG` you supply, e.g. one per thread.

```@example guide
leiden_clustering(g; seed=7).membership == leiden_clustering(g; seed=7).membership
```

```@example guide
rng = Xoshiro(7)
a = leiden_clustering(g; rng)
b = leiden_clustering(g; rng)   # continues the same stream, so it may differ from `a`
a.membership == leiden_clustering(g; seed=7).membership
```

Random streams are not guaranteed to be identical across Julia versions, so compare
modularity and the number of communities rather than exact memberships when results must
agree across machines.

Because Leiden is randomised, a single run is one sample. To keep the best of several:

```@example guide
best = argmax(r -> r.quality, [leiden_clustering(g; seed=s) for s in 1:10])
best.quality
```

## More Leiden passes

`n_iterations` (default 2, as in R/igraph) is the maximum number of full passes, each
started from the previous result. It stops early once a pass changes nothing. Extra passes
can improve quality at extra cost; a negative value runs until a pass changes nothing.

`initial_membership` starts the first pass from a given partition instead of singletons, for
example to refine the result of another method or of an earlier run. With `n_iterations=0`
it is returned unchanged (renumbered to `1:k`), with its quality:

```@example guide
start = leiden_clustering(g; seed=1).membership
leiden_clustering(g; initial_membership=start, n_iterations=-1, seed=2).quality
```

## Self-loops and requirements

- Graphs must be undirected; a directed graph throws an `ArgumentError`.
- A self-loop of weight `w` counts `2w` towards its vertex's weighted degree (igraph's
  convention), in both weighted and unweighted graphs.
- Multiple components are fine: communities never span components.
- A graph with no vertices returns an empty `Partition`; a graph with vertices but no
  edges returns one community per vertex, with quality `0.0`.
