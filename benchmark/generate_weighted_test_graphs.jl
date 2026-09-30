# generate_weighted_test_graphs.jl

using Random
using SparseArrays
using DataFrames
using Graphs
using SimpleWeightedGraphs
import CSV

# Helper: allow either a constant weight or a ()->Real sampler
_draw_weight(w) = w isa Function ? w() : w

"Erdős–Rényi G(n,p) weighted"
function weighted_er(n::Int, p::Float64; weight::Union{Real,Function}=1.0, seed::Integer=42)
    Random.seed!(seed)
    g = SimpleWeightedGraph(n)
    @inbounds for i in 1:n-1, j in (i+1):n
        if rand() < p
            add_edge!(g, i, j, _draw_weight(weight))
        end
    end
    return g
end

"2-block SBM with separate within/between probs and weights"
function weighted_sbm(n1::Int, n2::Int;
                      p_in::Float64=0.15, p_out::Float64=0.01,
                      w_in::Union{Real,Function}=1.0,
                      w_out::Union{Real,Function}=0.2,
                      seed::Integer=42)
    Random.seed!(seed)
    n = n1 + n2
    g = SimpleWeightedGraph(n)

    # block 1
    @inbounds for i in 1:n1-1, j in (i+1):n1
        if rand() < p_in
            add_edge!(g, i, j, _draw_weight(w_in))
        end
    end
    # block 2
    @inbounds for i in (n1+1):(n-1), j in (i+1):n
        if rand() < p_in
            add_edge!(g, i, j, _draw_weight(w_in))
        end
    end
    # between blocks
    @inbounds for i in 1:n1, j in (n1+1):n
        if rand() < p_out
            add_edge!(g, i, j, _draw_weight(w_out))
        end
    end
    return g
end

"Barabási–Albert preferential attachment (simple implementation), weighted"
function weighted_ba(n::Int, m::Int;
                     weight::Union{Real,Function}=1.0,
                     seed::Integer=42)
    @assert m ≥ 1 "m must be ≥ 1"
    @assert n ≥ m+1 "n must be at least m+1"

    Random.seed!(seed)
    g = SimpleWeightedGraph(n)

    # Start with a (m+1)-clique
    m0 = m + 1
    for i in 1:m0-1, j in (i+1):m0
        add_edge!(g, i, j, _draw_weight(weight))
    end

    degs = degree(g)
    totdeg = sum(degs)

    # Preferentially attach each new vertex t = m0+1..n with m edges
    for t in (m0+1):n
        targets = Int[]
        while length(targets) < m
            # sample a node proportional to degree
            r = rand() * totdeg
            acc = 0.0
            chosen = 1
            @inbounds for v in 1:(t-1) # only existing vertices
                acc += degs[v]
                if acc ≥ r
                    chosen = v
                    break
                end
            end
            if chosen != t && !(chosen in targets) # avoid self-loops & duplicates
                push!(targets, chosen)
            end
        end
        for v in targets
            add_edge!(g, t, v, _draw_weight(weight))
            degs[t] += 1
            degs[v] += 1
            totdeg += 2
        end
    end
    return g
end

"Watts–Strogatz small-world (ring lattice rewiring), weighted"
function weighted_ws(n::Int, k::Int, beta::Float64;
                     weight::Union{Real,Function}=1.0,
                     seed::Integer=42)
    @assert iseven(k) "k must be even"
    Random.seed!(seed)
    g = SimpleWeightedGraph(n)

    # start with k/2 neighbors on each side (ring lattice)
    halfk = k ÷ 2
    for i in 1:n
        for d in 1:halfk
            j = ((i - 1 + d) % n) + 1
            if !has_edge(g, i, j)
                add_edge!(g, i, j, _draw_weight(weight))
            end
        end
    end

    # rewire clockwise edges with prob beta
    for i in 1:n
        for d in 1:halfk
            j = ((i - 1 + d) % n) + 1
            if i < j && rand() < beta
                # remove current edge
                rem_edge!(g, i, j)
                # pick a new target avoiding self-loops and duplicates
                newj = i
                while newj == i || has_edge(g, i, newj)
                    newj = rand(1:n)
                end
                add_edge!(g, i, newj, _draw_weight(weight))
            end
        end
    end
    return g
end

"Axis-aligned 2D grid (4-neighborhood), weighted"
function weighted_grid_2d(rows::Int, cols::Int;
                          weight::Union{Real,Function}=1.0)
    n = rows * cols
    g = SimpleWeightedGraph(n)
    # map (r,c) -> 1..n
    idx(r,c) = (r-1)*cols + c

    for r in 1:rows, c in 1:cols
        u = idx(r,c)
        if r < rows
            v = idx(r+1,c)
            add_edge!(g, u, v, _draw_weight(weight))
        end
        if c < cols
            v = idx(r,c+1)
            add_edge!(g, u, v, _draw_weight(weight))
        end
    end
    return g
end

###############################################################################

graphs_w = Dict{String, SimpleWeightedGraph}()

graphs_w["w_er_1"]    = weighted_er(2_000, 0.0015)
graphs_w["w_sbm_1"]   = weighted_sbm(
    2_000, 2_000; p_in=0.02, p_out=0.002,
    w_in=1.0, w_out=0.2
)
graphs_w["w_ba_1"]    = weighted_ba(5_000, 3)
graphs_w["w_ws_1"]    = weighted_ws(2_000, 10, 0.1)
graphs_w["w_grid_2d"] = weighted_grid_2d(100, 100)

function save_edgelist(g::SimpleGraph, path::String)
    srcs = Int[]
    dsts = Int[]
    wts  = Int[]
    for e in edges(g)
        push!(srcs, src(e)); push!(dsts, dst(e)); push!(wts, 1)
    end
    CSV.write(path, DataFrame(src=srcs, dst=dsts, weight=wts))
end

graphs_w = Dict{String, Tuple{SimpleGraph,SparseMatrixCSC{Float64}}}()

# Save all
meta = DataFrame(name=String[], n=Int[], m=Int[])
for (name, g) in graphs_w
    path = "benchmark/data/$(name)_w.csv"
    save_edgelist(g, path)
    push!(meta, (name, nv(g), ne(g)))
end

CSV.write("benchmark/data/_meta_w.csv", meta)
println("Saved $(length(graphs_w)) graphs to benchmark/data/")
