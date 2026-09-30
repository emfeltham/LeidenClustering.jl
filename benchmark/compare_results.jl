# compare_results.jl

using CSV, DataFrames, Graphs
using LeidenClustering

function load_graph(path::String)
    df = CSV.read(path, DataFrame)
    n = maximum(vcat(df.src, df.dst))
    g = SimpleGraph(n)
    for r in eachrow(df); add_edge!(g, Int(r.src), Int(r.dst)); end
    return g
end

j = CSV.read("benchmark/julia/_summary_julia.csv", DataFrame)
r = CSV.read("benchmark/r/_summary_r.csv", DataFrame)

combined = innerjoin(j, r, on=[:name, :n, :m], makeunique=true)
combined.dQ = combined.modularity - combined.mod
rename!(
    combined, Dict(:mod => :Q_r, :modularity => :Q_julia,
    :communities => :k_julia, :communities_1 => :k_r,
    :julia_time_ms => :t_julia_ms, :r_time_ms => :t_r_ms)
)

combined.k_diff = combined.k_julia - combined.k_r;
combined.weighted = occursin.("w_", combined.name);

CSV.write("benchmark/_summary_combined.csv", combined)
println("Wrote combined summary to benchmark/_summary_combined.csv")

combined
