#!/bin/bash
# Download the frozen benchmark graphs (benchmark/data/) from the GitHub release
# `benchmark-data-v1` and verify their checksum.
#
# Usage: ./benchmark/fetch_data.sh [--force]
#   --force                   re-download even if benchmark/data/ already exists
#   BENCHMARK_DATA_TAR=path   use a local copy of the archive instead of downloading

set -euo pipefail

TAG="benchmark-data-v1"
ASSET="benchmark-data-v1.tar"
SHA256="ac22f8cd42b8150653bee7bc9ac81f1bf0690c4cebfde3b596689fe673977d0a"
URL="https://github.com/emfeltham/LeidenClustering.jl/releases/download/${TAG}/${ASSET}"

cd "$(dirname "$0")/.."   # repository root

if [ -f benchmark/data/_meta.csv ] && [ "${1:-}" != "--force" ]; then
    echo "benchmark/data/ already present (use --force to re-download)"
    exit 0
fi

archive="${BENCHMARK_DATA_TAR:-}"
if [ -z "$archive" ]; then
    archive="$(mktemp)"
    trap 'rm -f "$archive"' EXIT
    echo "Downloading $URL"
    curl -fL --retry 3 -o "$archive" "$URL"
fi

echo "$SHA256  $archive" | shasum -a 256 -c -
tar -xf "$archive" -C benchmark
echo "Extracted benchmark/data/ ($(ls benchmark/data | wc -l | tr -d ' ') files)"
