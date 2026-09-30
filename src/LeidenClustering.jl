# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Eric Martin Feltham

"""
    LeidenClustering

Community detection for undirected graphs: the Leiden and Louvain algorithms, after the
igraph C library, for any `Graphs.jl` graph (edge weights are read from
`SimpleWeightedGraph`).

# Exported
- [`leiden_clustering`](@ref): Leiden algorithm (well-connected communities), modularity or CPM
- [`louvain_clustering`](@ref): multilevel Louvain algorithm, modularity
- [`Partition`](@ref): the result type; [`ncommunities`](@ref) counts its communities

`LeidenClustering.modularity(g, membership; resolution=1.0)` is also available but not
exported, because `Graphs.modularity` has the same name.

# Example
```julia
using Graphs, LeidenClustering

g = smallgraph(:karate)
result = leiden_clustering(g; resolution=1.0, seed=42)
result.membership, result.quality
```
"""
module LeidenClustering

using Graphs
using Random
using SimpleWeightedGraphs
using SparseArrays

include("graph.jl")
include("workspace.jl")
include("partition.jl")
include("quality.jl")
include("result.jl")
include("louvain.jl")
include("leiden.jl")

export leiden_clustering, louvain_clustering, Partition, ncommunities

end # module
