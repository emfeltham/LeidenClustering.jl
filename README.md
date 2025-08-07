# LeidenClustering.jl
Leiden clustering on graphs.

**LeidenClustering.jl** is a pure Julia implementation of the [Leiden community detection algorithm](https://www.nature.com/articles/s41598-019-41695-z) for graphs.
It improves upon the Louvain method by ensuring that all communities are *well-connected*, avoiding the disconnected partitions that Louvain can produce.

The package supports weighted and unweighted graphs, arbitrary resolution parameters, and both the Louvain-like refinement and full Leiden refinement.

## Features

* Leiden algorithm with well-connectedness guarantee.
* Works with any `AbstractGraph` from [`Graphs.jl`](https://github.com/JuliaGraphs/Graphs.jl).
* Weighted or unweighted graphs.
* Adjustable resolution parameter for controlling community granularity.
* Deterministic: runs with `seed`.
* Modularity calculation and verification tools.
* Utilities for:

  * Displaying community structure
  * Sanity checks for partitions
  * Resolution sweeps

## Installation

The package is currently unregistered:

```julia
] add https://github.com/yourusername/LeidenClustering.jl
```

## Usage

```julia
using Graphs
using LeidenClustering

# Example graph: Zachary's Karate Club
g = Graphs.karateclub_graph()

# Run Leiden with full refinement
state = leiden(g; resolution=1.0, seed=42, use_refinement=true, use_true_leiden=true)

println("Communities found: ", length(state.partition.community_nodes))
println("Modularity: ", state.partition.modularity)
```

### Output partition

`state.partition` contains:

* `membership::Vector{Int}` — node → community mapping.
* `community_nodes::Dict{Int, Set{Int}}` — community → nodes.
* `modularity::Float64` — modularity of the current partition.
* `resolution::Float64` — resolution parameter used.

## API

### Core function

```julia
leiden(
    g::AbstractGraph;
    resolution::Float64 = 1.0,
    max_iterations::Int = 100,
    tol::Float64 = 1e-6,
    seed::Union{Int, Nothing} = nothing,
    use_refinement::Bool = true,
    use_true_leiden::Bool = false
) -> LeidenState
```

Arguments:

* `resolution` — resolution parameter γ. Lower → larger communities.
* `max_iterations` — maximum hierarchy levels.
* `tol` — modularity improvement tolerance for convergence.
* `seed` — RNG seed for reproducibility.
* `use_refinement` — whether to refine communities after local moving.
* `use_true_leiden` — if `true`, uses full Leiden refinement; if `false`, uses simple connectivity split.

### Utility functions

* `compute_modularity(state)` — compute modularity from current state.
* `verify_modularity(state)` — recompute modularity from scratch for validation.
* `show_communities(state)` — print summary of each community.
* `sanity_check(state)` — check for common consistency errors.

## Example: Resolution sweep

```julia
g = Graphs.karateclub_graph()

for res in [0.5, 1.0, 1.5]
    state = leiden(g; resolution=res, seed=42, use_refinement=true, use_true_leiden=true)
    println("Resolution: $res, Modularity: $(state.partition.modularity), Communities: $(length(state.partition.community_nodes))")
end
```

## Refinement

Note the refinement options on clustering.

* **Simple refinement** (`use_true_leiden=false`) only ensures communities are connected — effectively Louvain with a connectivity fix.
* **True refinement** (`use_true_leiden=true`) runs intra-community local moving to split communities into *well-connected* subcommunities before aggregation, following Traag et al. (2019).

## References

* V.A. Traag, L. Waltman, N.J. van Eck,
  *From Louvain to Leiden: guaranteeing well-connected communities*,
  *Scientific Reports* **9**, 5233 (2019).
  [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)
