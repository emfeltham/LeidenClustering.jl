# LeidenClustering.jl

[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://emfeltham.github.io/LeidenClustering.jl/dev/)
[![CI](https://github.com/emfeltham/LeidenClustering.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emfeltham/LeidenClustering.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/emfeltham/LeidenClustering.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emfeltham/LeidenClustering.jl)
[![License: GPL v3+](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](LICENSE)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Leiden and Louvain community detection for [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl)
graphs, in pure Julia. Both are derived from the [igraph](https://igraph.org) C implementations
(`leiden.c`, `louvain.c`) and support weighted graphs, a resolution parameter, and the modularity
and CPM objectives (Leiden).

- **Leiden** ([Traag et al. 2019](https://doi.org/10.1038/s41598-019-41695-z)) adds a refinement
  step to Louvain so communities are guaranteed to be well connected.
- **Louvain** ([Blondel et al. 2008](https://doi.org/10.1088/1742-5468/2008/10/P10008)) is the
  classic multilevel modularity optimiser. Communities may be internally disconnected.

## Installation

The package is not yet registered:

```julia
using Pkg
Pkg.add(url="https://github.com/emfeltham/LeidenClustering.jl")
```

Requires Julia 1.10 or later.

## Usage

```julia
using Graphs, LeidenClustering

g = smallgraph(:karate)

result = leiden_clustering(g; resolution=1.0, seed=42)
result.membership          # community id 1:k for every vertex
result.quality             # final modularity, ≈ 0.4198
ncommunities(result)      # 4

result = louvain_clustering(g; seed=42)

membership, quality = result   # a Partition also destructures like a tuple
```

Weighted graphs use their edge weights (or pass `weights=W`, any matrix with `W[i, j]` the
weight of edge `i–j`):

```julia
using SimpleWeightedGraphs
g = SimpleWeightedGraph(5)
add_edge!(g, 1, 2, 2.0); add_edge!(g, 2, 3, 1.0); add_edge!(g, 3, 4, 1.0); add_edge!(g, 4, 5, 2.0)
result = leiden_clustering(g)
result.membership   # [1, 1, 1, 2, 2]
```

## API

### `leiden_clustering(graph; resolution=1.0, beta=0.01, n_iterations=2, objective=:modularity, weights, node_weights, initial_membership, rng, seed)`

| Keyword | Meaning |
|---|---|
| `resolution` | γ. Higher gives more, smaller communities. |
| `beta` | Randomness of the refinement step (`0` = greedy). |
| `n_iterations` | Maximum number of full Leiden passes, each started from the previous result. Stops early once a pass changes nothing. Default 2, as in R/igraph; negative runs until nothing changes; `0` returns `initial_membership`. |
| `objective` | `:modularity`, or `:cpm` (Constant Potts Model). |
| `weights` | Edge weights, `weights[i, j]` for edge `i–j`, non-negative. Default `Graphs.weights(graph)`: the weights of a `SimpleWeightedGraph` (or any graph type defining it), else 1. |
| `node_weights` | CPM vertex weights (default 1 each). Not allowed with `:modularity`. |
| `initial_membership` | Starting partition, one positive id per vertex (default: singletons). |
| `rng` | Any `AbstractRNG`; defaults to the global RNG. |
| `seed` | Integer. Runs with a private `Xoshiro(seed)` instead of `rng`; the global RNG is not touched. |

### `louvain_clustering(graph; resolution=1.0, max_iterations=100, weights, rng, seed)`

Modularity only. `max_iterations` bounds the number of aggregation levels; `weights` as for Leiden.

### `Partition`

Both functions return a `Partition`:

| Field | Meaning |
|---|---|
| `membership` | Community id `1:k` of every vertex. |
| `quality` | Final quality of `membership` (modularity, or CPM for `objective=:cpm`). |
| `qualities` | Louvain: modularity after each level. Leiden: `[quality]`. |

`ncommunities(result)` gives `k`. Destructuring, `membership, quality = result`, works too.

### Other

- `LeidenClustering.modularity(graph, membership; resolution=1.0, weights)`: modularity of any partition
  of an undirected graph; community ids need not be consecutive. Not exported, because
  `Graphs.modularity` has the same name.

Both clustering functions require an undirected graph. A self-loop of weight `w` counts `2w`
toward its vertex's strength (igraph's convention).

## Documentation

A guide (weighted graphs, resolution, objectives, reproducibility), algorithm notes and the API
reference are published at <https://emfeltham.github.io/LeidenClustering.jl/dev/> (sources in `docs/`). Build them locally with

```bash
julia --project=docs docs/make.jl     # output in docs/build/, needs Julia 1.11+
```

## Validation

`benchmark/` compares Leiden with R/igraph's `cluster_leiden` (`resolution=1`, modularity) on
15 generated graphs (ER, BA, Watts–Strogatz, SBM, 2D grid, karate, weighted variants) with up to
10,000 vertices. Against R's saved results, the final modularity differs by 0.0011 on average and
0.0058 at most (Julia is higher on the largest gap); community counts are within a few. See
[`benchmark/README.md`](benchmark/README.md) to rerun. Leiden is randomised, so compare
modularity and community count, not memberships.

## Tests

```julia
using Pkg; Pkg.test()
```

The suite covers input validation, edge cases (empty and trivial graphs), determinism under a
seed, weighted graphs, self-loops, CPM, connectedness of Leiden communities, regression tests for
past bugs, and package hygiene via [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl). It is run
against Julia 1.10, 1.12 and 1.13.

## License and credits

LeidenClustering.jl is free software under the **GNU General Public License, version 3 or later**
(see [`LICENSE`](LICENSE)). It is derived from the C core of [igraph](https://igraph.org) (GPL-2.0-or-later),
by Gabor Csardi, Tamas Nepusz and the igraph development team; [`NOTICE.md`](NOTICE.md) lists which
igraph files each part comes from. A program that combines this package with your own code must itself be
distributed under GPL-compatible terms; using it to analyse data has no such requirement.

## References

- V. A. Traag, L. Waltman, N. J. van Eck, *From Louvain to Leiden: guaranteeing well-connected
  communities*, Scientific Reports **9**, 5233 (2019).
  [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)
- V. D. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre, *Fast unfolding of communities in
  large networks*, J. Stat. Mech. P10008 (2008).
  [doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)
- G. Csárdi, T. Nepusz, *The igraph software package for complex network research* (2006).
