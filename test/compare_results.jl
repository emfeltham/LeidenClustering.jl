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

j = CSV.read("test/julia/_summary_julia.csv", DataFrame)
r = CSV.read("test/r/_summary_r.csv", DataFrame)

combined = innerjoin(j, r, on=[:name, :n, :m], makeunique=true)
combined.dQ = combined.modularity - combined.mod
rename!(
    combined, Dict(:mod => :Q_r, :modularity => :Q_julia,
    :communities => :k_julia, :communities_1 => :k_r,
    :julia_time_ms => :t_julia_ms, :r_time_ms => :t_r_ms)
)

CSV.write("test/_summary_combined.csv", combined)
println("Wrote combined summary to data/_summary_combined.csv")
combined

# Optional sanity: compute Julia modularity using R’s partition
for row in eachrow(combined)
    name = row.name
    g = load_graph("data/graphs/$(name).csv")
    # Build a LeidenState with R partition
    state = LeidenClustering.initialize_state(g, resolution=1.0)
    rpart = CSV.read("data/r/$(name)_r_partition.csv", DataFrame)
    state.partition.membership = Int.(rpart.community)

    # rebuild structures and compute modularity in Julia
    LeidenClustering.rebuild_partition_structures!(state)
    Q_r_in_julia = state.partition.modularity

    println("[$(name)] Q_r_in_julia=$(round(Q_r_in_julia, digits=6)) vs Q_r_reported=$(round(row.Q_r, digits=6))")
end
