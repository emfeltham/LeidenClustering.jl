#!/bin/bash
# Run complete Julia vs R/igraph comparison pipeline
#
# Usage: ./benchmark/run_comparison.sh > ./benchmark/run_comparison.txt 2>&1

set -e  # Exit on error

echo "================================================================================"
echo "LEIDEN ALGORITHM: Julia vs R/igraph Comparison Pipeline"
echo "================================================================================"
echo ""

./benchmark/fetch_data.sh
echo ""

echo "Step 1/3: Running Julia Leiden benchmarks..."
echo "--------------------------------------------------------------------------------"
julia --project=benchmark benchmark/run_julia_leiden.jl
echo ""

echo "Step 2/3: Running R/igraph Leiden benchmarks..."
echo "--------------------------------------------------------------------------------"
Rscript benchmark/igraph_comparison.r
echo ""

echo "Step 3/3: Comparing results..."
echo "--------------------------------------------------------------------------------"
julia --project=benchmark benchmark/compare_results.jl
echo ""

echo "================================================================================"
echo "Displaying comparison results..."
echo "================================================================================"
echo ""

# Display formatted comparison
julia --project=benchmark -e '
using CSV, DataFrames, Printf, Statistics

df = CSV.read("benchmark/_summary_combined.csv", DataFrame)

println("="^80)
println("LEIDEN ALGORITHM: Julia vs R/igraph Comparison")
println("="^80)

println("\nKey Metrics:")
println("  - dQ: Modularity difference (Julia - R)")
println("  - k_diff: Community count difference (Julia - R)")
println("  - Negative dQ means R found higher modularity")
println("\n" * "-"^80)

for row in eachrow(df)
    @printf("%-12s | n=%5d m=%6d | Julia: k=%3d Q=%.4f | R: k=%3d Q=%.4f | ΔQ=%+.4f Δk=%+4d\n",
        row.name, row.n, row.m,
        row.k_julia, row.Q_julia,
        row.k_r, row.Q_r,
        row.dQ, row.k_diff)
end

println("-"^80)
println("\nSummary Statistics:")
@printf("  Average |ΔQ|: %.4f\n", mean(abs.(df.dQ)))
@printf("  Max |ΔQ|: %.4f (%s)\n", maximum(abs.(df.dQ)), df.name[argmax(abs.(df.dQ))])
@printf("  Graphs with similar modularity (|ΔQ| < 0.01): %d/%d\n",
    count(abs.(df.dQ) .< 0.01), nrow(df))
@printf("  Graphs with exact same k: %d/%d\n", count(df.k_diff .== 0), nrow(df))
@printf("  Average |Δk|: %.1f communities\n", mean(abs.(df.k_diff)))

println("\n" * "="^80)
println("CONCLUSION:")
println("="^80)
avg_dq = mean(abs.(df.dQ))
if avg_dq < 0.01
    println("  ✓✓ Julia implementation has EXCELLENT agreement with R/igraph")
    println("  The implementations produce nearly identical results.")
elseif avg_dq < 0.05
    println("  ✓ Julia implementation has GOOD agreement with R/igraph")
    println("  Small differences are expected due to randomization and implementation details.")
else
    println("  ⚠ Julia implementation differs significantly from R/igraph")
end
@printf("\n  Average modularity difference: %.4f\n", avg_dq)
println("\n  Results saved to: benchmark/_summary_combined.csv")
println("="^80 * "\n")
'

echo "✓ Comparison complete!"
echo ""
