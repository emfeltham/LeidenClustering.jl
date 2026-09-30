# Notice

LeidenClustering.jl is Copyright (C) 2025 Eric Martin Feltham and is licensed under the
GNU General Public License, version 3 or (at your option) any later version
(`GPL-3.0-or-later`). See [`LICENSE`](LICENSE).

## Derived from igraph

This package began as a Julia translation of, and remains derived from, parts of the C core of the
[igraph](https://igraph.org) network analysis library, which is licensed under the GNU General Public
License, version 2 or (at your option) any later version. That is why this package is distributed under
the GPL and not a permissive licence. The following igraph source files were translated or consulted:

| igraph file | igraph copyright notice | Julia code derived from it |
|---|---|---|
| `src/community/leiden.c` | Copyright (C) 2007-2012 Gabor Csardi <csardi.gabor@gmail.com> | `src/leiden.jl` |
| `src/community/louvain.c` | Copyright (C) 2007-2020 The igraph development team | `src/louvain.jl` |
| `src/properties/degrees.c` | Copyright (C) 2005-2023 The igraph development team | weighted-degree and self-loop handling in `src/graph.jl`, `src/workspace.jl` and `src/louvain.jl` |
| `src/community/modularity.c` | Copyright (C) 2007-2020 The igraph development team | `src/quality.jl` |

We are grateful to Gabor Csardi, Tamas Nepusz and the igraph development team, and to everyone who has
contributed to igraph, for the reference implementations this package builds on. The package has since been
substantially restructured and extended (a shared sparse-matrix core, a queue-based Louvain local-moving
step, a result type, explicit random number generators), and it is not endorsed by or affiliated with the
igraph project. Where it departs from igraph, see the "Differences from igraph" section of the
documentation.

If you use igraph itself, please cite: G. Csardi and T. Nepusz, *The igraph software package for complex
network research*, InterJournal Complex Systems, 1695 (2006).

## Algorithms

The algorithms are due to their authors; please cite them when you use the results:

- V. A. Traag, L. Waltman, N. J. van Eck, *From Louvain to Leiden: guaranteeing well-connected
  communities*, Scientific Reports **9**, 5233 (2019). <https://doi.org/10.1038/s41598-019-41695-z>
- V. D. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre, *Fast unfolding of communities in large
  networks*, Journal of Statistical Mechanics: Theory and Experiment P10008 (2008).
  <https://doi.org/10.1088/1742-5468/2008/10/P10008>
- M. E. J. Newman, *Modularity and community structure in networks*, PNAS **103**, 8577 (2006).

## Validation

The scripts in `benchmark/` compare results with R's `igraph` package (also GPL). No code or data from it
is included in this repository.
