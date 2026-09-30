# Validation

The `benchmark/` directory compares Leiden with R/igraph's `cluster_leiden`
(`resolution=1`, modularity objective) on 15 frozen graphs: Erdős–Rényi, Barabási–Albert,
Watts–Strogatz, stochastic block model, a 2D grid and Zachary's karate club, plus weighted
variants, with up to 10 000 vertices.

Against R's saved results, the final modularity differs by **0.0011 on average and 0.0058 at
most**, and community counts are within a few. Because both methods are randomised the two
values are not expected to be identical; what matters is that neither is systematically
better.

The test suite (`test/runtests.jl`) additionally checks that

- the modularity a function reports equals the modularity recomputed from its membership;
- Leiden's communities are connected;
- planted communities are recovered, with and without self-loops and edge weights;
- results are deterministic under a `seed` and never touch the global random number generator;
- the package passes [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl)'s hygiene checks on
  Julia 1.10, 1.12 and 1.13.

The frozen benchmark graphs are published as the `benchmark-data-v1` release asset, not stored in the
repository; `benchmark/fetch_data.sh` downloads and verifies them. To re-run the comparison with R, see
`benchmark/README.md`. It needs R with the `igraph` and
`readr` packages and uses its own Julia environment (`benchmark/Project.toml`).
