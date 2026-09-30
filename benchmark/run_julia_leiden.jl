# run_julia_leiden.jl

using Random
using Graphs, SimpleWeightedGraphs, CSV, DataFrames, Printf
using LeidenClustering

const RES = 1.0
const USE_TRUE_LEIDEN = true

function load_graph(path::String)
    df = CSV.read(path, DataFrame)
    n = maximum(vcat(df.src, df.dst))
    if :weight in propertynames(df)
        g = SimpleWeightedGraph(n)
        for row in eachrow(df)
            add_edge!(g, Int(row.src), Int(row.dst), Float64(row.weight))
        end
        return g
    end
    g = SimpleGraph(n)
    for row in eachrow(df)
        add_edge!(g, Int(row.src), Int(row.dst))
    end
    return g
end

# Best of `reps` runs after one warm-up call (seconds)
function best_time(f; reps=3)
    f()
    return minimum(begin t0 = time_ns(); f(); (time_ns() - t0) / 1e9 end for _ in 1:reps)
end

isfile("benchmark/data/_meta.csv") ||
    error("benchmark/data/ is missing: run ./benchmark/fetch_data.sh from the package root")

meta = CSV.read("benchmark/data/_meta.csv", DataFrame);
meta_w = CSV.read("benchmark/data/_meta_w.csv", DataFrame);
meta = vcat(meta, meta_w);

# initialize storage
summ = DataFrame(
    name=String[], n=Int[], m=Int[], julia_time_ms=Float64[],
    communities=Int[], modularity=Float64[]
);

Random.seed!(08540)

# row = collect(eachrow(meta))[1];
for row in eachrow(meta)
    name = row.name
    g = load_graph("benchmark/data/$(name).csv.gz")

    t = best_time(() -> leiden_clustering(g; resolution=RES, objective=:modularity, seed=1))

    result = leiden_clustering(g; resolution=RES, objective=:modularity)
    n_clusts = ncommunities(result)
    quality = result.quality

    push!(
        summ, (name, nv(g), ne(g), t*1000, n_clusts, quality)
    )
end

s1 = deepcopy(summ)

mkpath("benchmark/julia")
CSV.write("benchmark/julia/_summary_julia.csv", summ)
