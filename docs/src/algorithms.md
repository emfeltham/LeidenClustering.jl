# Algorithms

## Modularity and CPM

Both algorithms maximise a quality function of the form

```math
Q = \frac{1}{2m} \sum_{ij} \left( A_{ij} - \gamma \, n_i n_j \right) \delta(\sigma_i, \sigma_j),
```

where ``A`` is the weighted adjacency matrix, ``\sigma_i`` the community of vertex ``i``, and
``\gamma`` the `resolution`. The two objectives differ in the vertex weights ``n_i``:

- **Modularity**: ``n_i = k_i / 2m`` with ``k_i`` the weighted degree, so ``\gamma n_i n_j``
  is the expected weight of an edge between ``i`` and ``j`` in a random graph with the same
  degrees.
- **CPM**: ``n_i = 1``, so ``\gamma`` is a density threshold that does not depend on the rest
  of the graph.

## Louvain

[`louvain_clustering`](@ref) repeats two steps until modularity stops improving:

1. **Local moving.** Visit vertices in random order and move each to the neighbouring
   community with the largest positive gain in modularity.
2. **Aggregation.** Contract every community into one vertex. Weight inside a community
   becomes a self-loop, so modularity is unchanged by contraction.

## Leiden

[`leiden_clustering`](@ref) adds a refinement step so communities are always connected:

1. **Fast local moving.** As in Louvain, but with a work queue: a vertex is revisited only
   after a neighbour moved.
2. **Refinement.** Inside each community, start from singletons and merge vertices that are
   well connected to the community, choosing the target at random with probability
   proportional to ``\exp(\Delta Q / \beta)``. `beta` controls the randomness; `beta=0` is
   greedy.
3. **Aggregation.** Contract the *refined* communities. The contracted vertices start out in
   the communities from step 1, so the coarse graph can still merge them.

Passes repeat on the aggregated graph until every community is a single vertex, and the whole
procedure is repeated up to `n_iterations` times from the previous result (the first run
starts from `initial_membership`, or from singletons).

## Implementation

The graph is held as a symmetric sparse matrix (CSC) with a separate vector of self-loop
weights. Both algorithms share a sparse accumulator for the weight a vertex has to each
neighbouring community, a counting sort to group vertices by community, and a contraction
routine that builds the aggregated CSC directly without sorting triplets. Everything derived
from the original graph (the total weight ``2m``) is computed once, because Leiden's
aggregation drops the weight inside communities.

## Differences from igraph

The package started as a translation of igraph's `leiden.c` and `louvain.c` and follows them
closely, with these deliberate exceptions:

- **Louvain uses a work queue for local moving**, not igraph's repeated full sweeps over all
  vertices. Both stop at a local optimum under the same "strictly positive gain" rule, but the
  queue visits far fewer vertices: 4–34× faster on the benchmark graphs and 16× on a
  100 000-vertex graph. Modularity is within about 0.002 of the sweep version.
- **Random streams differ.** Julia's random number generators, and the order in which
  candidate communities are visited, differ from igraph's, so memberships are not
  reproducible against R or C for the same seed. Compare modularity and the number of
  communities instead.
- **Results are a [`Partition`](@ref)** rather than igraph's membership vector.
- **Directed graphs are rejected.** igraph's Leiden also requires undirected graphs.

## References

- V. A. Traag, L. Waltman, N. J. van Eck, *From Louvain to Leiden: guaranteeing well-connected
  communities*, Scientific Reports **9**, 5233 (2019).
  [doi:10.1038/s41598-019-41695-z](https://doi.org/10.1038/s41598-019-41695-z)
- V. D. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre, *Fast unfolding of communities in
  large networks*, J. Stat. Mech. P10008 (2008).
  [doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)
- M. E. J. Newman, *Modularity and community structure in networks*, PNAS **103**, 8577 (2006).
- G. Csárdi, T. Nepusz, *The igraph software package for complex network research*,
  InterJournal Complex Systems 1695 (2006).
