# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Eric Martin Feltham
# Self-loop and strength conventions follow igraph (degrees.c, Copyright (C) 2005-2023 The igraph
# development team, GPL-2.0-or-later). See NOTICE.md.

# Internal weighted-graph representation shared by Louvain and Leiden.

"""
    Adjacency

Undirected weighted graph as a symmetric, zero-diagonal sparse matrix `A` (CSC) plus a
separate vector of self-loop weights. A loop of weight `w` contributes `2w` to its
vertex's strength (igraph's convention). Aggregation keeps intra-community weight as
loops, so the same type serves every level of the multilevel algorithms.
"""
struct Adjacency
    A::SparseMatrixCSC{Float64,Int}
    loops::Vector{Float64}
end

Graphs.nv(adj::Adjacency) = size(adj.A, 1)

"""
    Adjacency(graph::AbstractGraph, weights=Graphs.weights(graph))

Build from any undirected Graphs.jl graph. `weights[i, j]` is read for every edge, so any
`AbstractMatrix` works (for a loop, `weights[i, i]` is its weight `w`). The default,
`Graphs.weights(graph)`, gives edge weights for `SimpleWeightedGraph` and other weighted
graph types, and 1 for every edge otherwise. Throws if a weight is negative or not finite.
"""
function Adjacency(graph::AbstractGraph, weights=Graphs.weights(graph))
    adj = if weights isa Graphs.DefaultDistance && !(graph isa AbstractSimpleWeightedGraph)
        # adjacency_matrix stores 2 on the diagonal for a loop
        split_diagonal(SparseMatrixCSC{Float64,Int}(adjacency_matrix(graph)), 0.5)
    elseif graph isa AbstractSimpleWeightedGraph && weights === graph.weights
        # weights[i, i] holds the loop weight w itself
        split_diagonal(SparseMatrixCSC{Float64,Int}(graph.weights), 1.0)
    else
        size(weights) == (nv(graph), nv(graph)) ||
            throw(ArgumentError("weights must be an nv(graph) × nv(graph) matrix"))
        edge_weight_matrix(graph, weights)
    end
    check_weights(adj)
    return adj
end

# Symmetric weight matrix of the edges of `graph`, loops on the diagonal (weight w).
function edge_weight_matrix(graph::AbstractGraph, weights)
    m = ne(graph)
    I = Vector{Int}(undef, 2m)
    J = Vector{Int}(undef, 2m)
    V = Vector{Float64}(undef, 2m)
    k = 0
    for e in edges(graph)
        i, j = src(e), dst(e)
        w = Float64(weights[i, j])
        k += 1
        I[k], J[k], V[k] = i, j, w
        if i != j
            k += 1
            I[k], J[k], V[k] = j, i, w
        end
    end
    resize!(I, k); resize!(J, k); resize!(V, k)
    n = nv(graph)
    return split_diagonal(sparse(I, J, V, n, n), 1.0)
end

function check_weights(adj::Adjacency)
    ok(w) = isfinite(w) && w >= 0
    (all(ok, nonzeros(adj.A)) && all(ok, adj.loops)) ||
        throw(ArgumentError("Edge weights must be non-negative and finite"))
    return nothing
end

# Move the diagonal of `M` (scaled by `loopscale`) into `loops`.
function split_diagonal(M::SparseMatrixCSC{Float64,Int}, loopscale::Float64)
    n = size(M, 1)
    rows = rowvals(M)
    vals = nonzeros(M)
    loops = zeros(Float64, n)
    colptr = Vector{Int}(undef, n + 1)
    rowval = Vector{Int}(undef, nnz(M))
    nzval = Vector{Float64}(undef, nnz(M))
    pos = 1
    colptr[1] = 1
    @inbounds for j in 1:n
        for p in nzrange(M, j)
            i = rows[p]
            if i == j
                loops[j] += loopscale * vals[p]
            else
                rowval[pos] = i
                nzval[pos] = vals[p]
                pos += 1
            end
        end
        colptr[j + 1] = pos
    end
    resize!(rowval, pos - 1)
    resize!(nzval, pos - 1)
    return Adjacency(SparseMatrixCSC(n, n, colptr, rowval, nzval), loops)
end

"""
    strengths(adj) -> Vector{Float64}

Weighted degree of every vertex, self-loops counted twice.
"""
function strengths(adj::Adjacency)
    A = adj.A
    vals = nonzeros(A)
    s = Vector{Float64}(undef, size(A, 1))
    @inbounds for j in eachindex(s)
        t = 2 * adj.loops[j]
        for p in nzrange(A, j)
            t += vals[p]
        end
        s[j] = t
    end
    return s
end
