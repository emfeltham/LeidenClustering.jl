# Benchmarks

Validation of LeidenClustering.jl against R/igraph. Not part of the test suite (`test/runtests.jl`).
Run everything from the package root.

## Layout

| Path | Purpose |
|---|---|
| `fetch_data.sh` | Downloads and checksums the frozen inputs into `data/` from the `benchmark-data-v1` GitHub release. |
| `data/` | Frozen benchmark inputs (gzipped edge lists `<name>.csv.gz` and `w_<name>.csv.gz`, and the small `_meta*.csv` index files). Not in git: fetched by `fetch_data.sh`. |
| `run_julia_leiden.jl` | Runs `leiden_clustering` on every graph in `data/`; writes `julia/_summary_julia.csv`. |
| `igraph_comparison.r` | Runs R/igraph `cluster_leiden` on the same graphs; writes `r/_summary_r.csv` and `data/r/`. |
| `compare_results.jl` | Joins the two summaries into `_summary_combined.csv`. |
| `run_comparison.sh` | Runs all three steps and prints a report. |
| `generate_*_graphs.jl` | Generators used to create `data/`. |

`data/`, `julia/`, `r/`, `_summary_combined.csv` and `run_comparison.txt` are generated and gitignored.

## Usage

The benchmarks use their own environment (`benchmark/Project.toml`, Julia ≥ 1.11 for `[sources]`).
First time: `./benchmark/fetch_data.sh` (about 1.5 MB, checksum-verified) and
`julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'`. `run_comparison.sh` fetches the data itself.

```bash
./benchmark/run_comparison.sh                 # full pipeline (needs R with igraph, readr)
julia --project=benchmark benchmark/run_julia_leiden.jl   # Julia side only
```

## Notes

- The release asset is the source of truth (about 1.5 MB gzipped; `CSV.jl` and `readr` read `.csv.gz` directly). The generators do not reproduce it exactly (`ws_1`, `ws_2` and
  the `_meta*` files differ under current Graphs.jl, and the weighted generator is not working),
  so do not regenerate over it.
- The generators write plain `.csv`; gzip their output (`gzip -9 -n`) to match `data/`. A changed data set needs a new release asset
  (`benchmark-data-v2`) and a new checksum in `fetch_data.sh`; never replace an existing asset.
- Graph edge lists carry a `weight` column; `w_*` graphs are loaded as weighted graphs.
- Leiden is randomised: compare modularity and community count, not memberships.
