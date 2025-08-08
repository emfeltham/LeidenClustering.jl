# generate_test_graphs.jl

using Graphs, Random, CSV, DataFrames, Printf
using LeidenClustering

function simple_sbm(n1::Int, n2::Int; p_in::Float64=0.15, p_out::Float64=0.01, seed=42)
    Random.seed!(seed)
    n = n1 + n2
    g = SimpleGraph(n)
    for i in 1:n1-1, j in (i+1):n1
        rand() < p_in && add_edge!(g, i, j)
    end
    for i in n1+1:n-1, j in (i+1):n
        rand() < p_in && add_edge!(g, i, j)
    end
    for i in 1:n1, j in (n1+1):n
        rand() < p_out && add_edge!(g, i, j)
    end
    return g
end

function save_edgelist(g::SimpleGraph, path::String)
    srcs = Int[]
    dsts = Int[]
    wts  = Int[]
    for e in edges(g)
        push!(srcs, src(e)); push!(dsts, dst(e)); push!(wts, 1)
    end
    CSV.write(path, DataFrame(src=srcs, dst=dsts, weight=wts))
end

graphs = Dict{String, SimpleGraph}()

# 10 graphs (mix of types)

Random.seed!(08540)

graphs["karate"]         = LeidenClustering.karateclub_graph()
graphs["er_1"]           = erdos_renyi(2_000, 0.0015) # ~3k edges
graphs["er_2"]           = erdos_renyi(5_000, 0.0009)
graphs["ba_1"]           = barabasi_albert(5_000, 3)
graphs["ba_2"]           = barabasi_albert(10_000, 5)
graphs["ws_1"]           = watts_strogatz(2_000, 10, 0.1) # k=10
graphs["ws_2"]           = watts_strogatz(5_000, 12, 0.05)
graphs["sbm_1"]          = simple_sbm(2_000, 2_000; p_in=0.02, p_out=0.002)
graphs["sbm_2"]          = simple_sbm(3_000, 3_000; p_in=0.015, p_out=0.0015)
graphs["grid_2d"]        = grid([100, 100]) |> SimpleGraph # 10k nodes


graphs_w = Dict{String, Tuple{SimpleGraph,SparseMatrixCSC{Float64}}}()

graphs_w["w_er_1"]   = weighted_er(2_000, 0.0015)
graphs_w["w_sbm_1"]  = weighted_sbm(2_000, 2_000; p_in=0.02, p_out=0.002,
                                    w_in=1.0, w_out=0.2)
graphs_w["w_ba_1"]   = weighted_ba(5_000, 3)
graphs_w["w_ws_1"]   = weighted_ws(2_000, 10, 0.1)
graphs_w["w_grid_2d"] = weighted_grid_2d(100, 100)

# Save all
meta = DataFrame(name=String[], n=Int[], m=Int[])
for (name, g) in graphs
    path = "data/$(name).csv"
    save_edgelist(g, path)
    push!(meta, (name, nv(g), ne(g)))
end

CSV.write("data/_meta.csv", meta)
println("Saved $(length(graphs)) graphs to data/graphs/")
