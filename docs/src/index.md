# LeidenClustering.jl

Community detection for [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) graphs, in pure
Julia, with two algorithms:

- **Leiden** ([`leiden_clustering`](@ref)) ([Traag et al. 2019](https://doi.org/10.1038/s41598-019-41695-z))
  refines Louvain's communities so that they are guaranteed to be connected. It supports the
  modularity and Constant Potts Model (CPM) objectives.
- **Louvain** ([`louvain_clustering`](@ref)) ([Blondel et al. 2008](https://doi.org/10.1088/1742-5468/2008/10/P10008))
  is the classic multilevel modularity optimiser, and is the faster of the two.

Both work on unweighted and weighted undirected graphs, take a resolution parameter, and are
derived from the [igraph](https://igraph.org) C implementations. Leiden's modularity agrees with
R/igraph's `cluster_leiden` to within 0.006 on every graph we compared (see
[Validation](validation.md)).

## License

LeidenClustering.jl is licensed under the **GNU General Public License, version 3 or later**. It is
derived from the C core of [igraph](https://igraph.org), which is GPL-2.0-or-later, and credits
Gabor Csardi, Tamas Nepusz and the igraph development team; see `NOTICE.md` in the repository for the
file-by-file attribution. A program that combines the package with other code must itself be distributed
under GPL-compatible terms.

## Installation

The package is not yet registered:

```julia
using Pkg
Pkg.add(url="https://github.com/emfeltham/LeidenClustering.jl")
```

It requires Julia 1.10 or later.

## Quick start

```@example quick
using Graphs, LeidenClustering

g = smallgraph(:karate)
result = leiden_clustering(g; seed=42)
```

The result is a [`Partition`](@ref):

```@example quick
result.membership
```

```@example quick
(k = ncommunities(result), modularity = result.quality)
```

Continue with the [Guide](guide.md) for weighted graphs, resolution, objectives and
reproducibility, or [Algorithms](algorithms.md) for how the methods work and where this
package differs from igraph.
