# Maintainer notes

Internal reminders; not user documentation (that lives in `docs/` and `README.md`).

## CI, docs and badges

Workflows in `.github/workflows/` (from the package template): `CI.yml` (tests on Julia 1.10, the current release
and prerelease on Linux; the release also on macOS/aarch64 and Windows; prerelease may fail; coverage uploaded to
Codecov), `Documentation.yml` (builds the docs on every push/PR, deploys `dev` from `main` and versioned docs from
tags), `TagBot.yml` and `CompatHelper.yml`. The README badges point at them. Do not overwrite these files
without looking at them first (`git diff .github`).

One-time setup on GitHub; the workflows have not run yet, so the badges show "no status" until they do:

- [ ] **Pages:** Settings -> Pages -> deploy from branch `gh-pages` (created by the first docs run on `main`).
- [ ] **Codecov:** enable the repo at codecov.io and add its upload token as the Actions secret `CODECOV_TOKEN`
      (`fail_ci_if_error: false`, so a missing token only leaves the coverage badge empty).
- [ ] **Docs deploy key (optional):** `Documentation.yml` uses `GITHUB_TOKEN`, and also `DOCUMENTER_KEY` if set.
- [ ] **Docs badge:** it links to `dev`. After the first `vX.Y.Z` release add a `stable` badge:
      `https://img.shields.io/badge/docs-stable-blue.svg` -> `https://emfeltham.github.io/LeidenClustering.jl/stable/`.
- [ ] **TagBot/CompatHelper** are only useful once the package is registered in the General registry.
- The docs job uses Julia `1` because `docs/Project.toml` relies on `[sources]` (needs >= 1.11).

## Licensing

The package is **GPL-3.0-or-later** (`LICENSE`), not MIT, because it is derived from igraph's C core
(GPL-2.0-or-later, "or later" is what makes GPL-3 possible). Attribution lives in `NOTICE.md` and in the
`SPDX`/"Derived from" headers at the top of each `src/` file.

- Keep the header on every new `src/` file; add the igraph source to `NOTICE.md` if new code is derived from
  another igraph file.
- The `igraph/` directory holds GPL C sources for reference only. It is gitignored; never commit it.
- Every commit before the relicensing (`main` at `e5716e2` and earlier) carries an MIT `LICENSE`. The code they contain was already derived from igraph, so the GPL terms are the ones that apply
  to it; do not tag or register any release from those commits.
- Before registering: confirm `LICENSE` is detected as GPL-3.0 on GitHub, and that the docs page and README
  state the licence.

## Benchmark data is a release asset, not in git

The frozen benchmark graphs (`benchmark/data/`, 17 files, ~1.5 MB) are **not committed**. They are
attached to the GitHub release `benchmark-data-v1` as `benchmark-data-v1.tar`, and
`benchmark/fetch_data.sh` downloads them and checks a pinned SHA-256.

Why: they must stay frozen (the generators in `benchmark/generate_*.jl` do **not** reproduce them under
current Graphs.jl, so they cannot be regenerated), but Pkg downloads the whole repo tree for every
user, so they should not ship with the package.

| Thing | Where |
|---|---|
| Release / tag | `benchmark-data-v1` on `emfeltham/LeidenClustering.jl` (created with `--latest=false`) |
| Asset | `benchmark-data-v1.tar` |
| Pinned checksum | `SHA256=` in `benchmark/fetch_data.sh` (currently `ac22f8cd…977d0a`) |
| Git rule | `/benchmark/data/` and `*.tar` in `.gitignore` |
| Local copy of the archive | `benchmark-data-v1.tar` in the repo root (ignored). **Keep a copy outside the repo too.** |

### Use it

```bash
./benchmark/fetch_data.sh                          # skips if benchmark/data/ exists; --force to redo
BENCHMARK_DATA_TAR=/path/to/archive.tar ./benchmark/fetch_data.sh   # no network
./benchmark/run_comparison.sh                      # fetches first, then runs Julia vs R
```

### Rules

- **Never replace or edit an existing asset.** The checksum is pinned in the script, so replacing it
  breaks every checkout of every older commit. A changed data set becomes `benchmark-data-v2`, with its
  own release and a new `TAG`, `ASSET` and `SHA256` in `fetch_data.sh`.
- The edge lists are gzipped (`gzip -9 -n`) and read as `.csv.gz` by CSV.jl and `readr`. The small
  `_meta*.csv` index files stay plain.
- If the release is deleted, `fetch_data.sh` fails with a 404. Recreate it from the local archive
  (below); the checksum will still match.

### Build the archive

Run from the repo root, with `benchmark/data/` populated. The flags matter: without them macOS adds
`.DS_Store` and `._*` files, which change the checksum.

```bash
COPYFILE_DISABLE=1 tar --exclude='.DS_Store' --exclude='._*' --exclude='data/r' \
    -cf benchmark-data-v1.tar -C benchmark data
tar -tf benchmark-data-v1.tar          # expect data/, data/_meta.csv, data/_meta_w.csv and 15 *.csv.gz
shasum -a 256 benchmark-data-v1.tar    # must equal SHA256 in benchmark/fetch_data.sh
```

Tar files are not byte-reproducible across machines (file order and metadata vary). If you rebuild
the archive rather than reuse the original, the checksum will probably differ; either reuse the
original file or publish it as a new version.

### Publish

```bash
gh release create benchmark-data-v1 benchmark-data-v1.tar \
  --repo emfeltham/LeidenClustering.jl --target main --latest=false \
  --title "Benchmark data v1" \
  --notes "Frozen benchmark graphs for benchmark/. Fetched by benchmark/fetch_data.sh; not part of the package."
```

Then check from a clean clone: `./benchmark/fetch_data.sh --force` should print `OK` and extract 17 files.

### Status

- [ ] Release `benchmark-data-v1` created (as of writing it had **not** been published; until it is,
      use `BENCHMARK_DATA_TAR`)
- [ ] Copy of `benchmark-data-v1.tar` stored outside the repo
- [ ] R side of the pipeline (`benchmark/igraph_comparison.r`, reads `.csv.gz`) run once against the
      release data (no `Rscript` was available when it was changed)
