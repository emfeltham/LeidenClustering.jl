# run_julia_leiden.jl

using Revise
using BenchmarkTools
using Random
using Graphs, CSV, DataFrames, Printf
using LeidenClustering

const RES = 1.0
const USE_TRUE_LEIDEN = true

function load_graph(path::String)
    df = CSV.read(path, DataFrame)
    n = maximum(vcat(df.src, df.dst))
    g = SimpleGraph(n)
    for row in eachrow(df)
        add_edge!(g, Int(row.src), Int(row.dst))
    end
    return g
end

meta = CSV.read("data/_meta.csv", DataFrame)

# initialize storage
summ = DataFrame(
    name=String[], n=Int[], m=Int[], julia_time_ms=Float64[],
    communities=Int[], modularity=Float64[]
);

Random.seed!(08540)

# row = meta[1,:]
for row in eachrow(meta)
    name = row.name
    g = load_graph("data/$(name).csv")

    t = @belapsed leiden(
        $g; resolution=$RES,
        use_refinement=true,
        use_true_leiden=$USE_TRUE_LEIDEN
    )
    state = leiden(
        g; resolution=RES,
        use_refinement=true,
        use_true_leiden=true
    );

    # CSV.write(
    #     "test/julia/$(name)_julia_partition.csv",
    #     DataFrame(node=1:nv(g), community=state.partition.membership)
    # )

    push!(
        summ, (name, nv(g), ne(g), t*1000, length(state.partition.community_nodes),
        state.partition.modularity)
    )
    println("Julia done: $(name) | Q=$(round(state.partition.modularity, digits=6)) | k=$(length(state.partition.community_nodes)) | $(round(t*1000, digits=1)) ms")
end

summ

CSV.write("test/julia/_summary_julia.csv", summ)
println("Wrote Julia summary to data/julia/_summary_julia.csv")
