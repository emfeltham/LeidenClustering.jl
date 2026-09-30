# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Eric Martin Feltham

# Operations on partitions: renumbering, grouping and contraction.

"""
    renumber!(membership) -> k

Renumber communities to `1:k` in order of first appearance and return `k`.
"""
function renumber!(membership::Vector{Int})
    isempty(membership) && return 0
    map = zeros(Int, maximum(membership))
    k = 0
    @inbounds for i in eachindex(membership)
        c = membership[i]
        if map[c] == 0
            k += 1
            map[c] = k
        end
        membership[i] = map[c]
    end
    return k
end

"""
    group_members(membership, k) -> (start, members)

Counting sort of the vertices by community (`membership` uses ids `1:k`): the members of
community `c` are `members[start[c]:start[c+1]-1]`.
"""
function group_members(membership::Vector{Int}, k::Int)
    start = zeros(Int, k + 1)
    @inbounds for c in membership
        start[c + 1] += 1
    end
    start[1] = 1
    @inbounds for c in 1:k
        start[c + 1] += start[c]
    end
    next = start[1:k]
    members = Vector{Int}(undef, length(membership))
    @inbounds for v in eachindex(membership)
        c = membership[v]
        members[next[c]] = v
        next[c] += 1
    end
    return start, members
end

"""
    contract(A, membership, k) -> (B, internal)

Contract each community of `membership` (ids `1:k`) into one vertex. `B` is the `k×k`
symmetric zero-diagonal weight matrix between communities; `internal[c]` is the total
weight of the edges inside community `c` (each edge counted once). Builds the CSC
structure directly, one column per community.
"""
function contract(A::SparseMatrixCSC{Float64,Int}, membership::Vector{Int}, k::Int)
    rows = rowvals(A)
    vals = nonzeros(A)
    start, members = group_members(membership, k)
    acc = NeighborAccumulator(k)
    colptr = Vector{Int}(undef, k + 1)
    rowval = Vector{Int}(undef, nnz(A))
    nzval = Vector{Float64}(undef, nnz(A))
    internal = zeros(Float64, k)
    pos = 1
    colptr[1] = 1
    @inbounds for c in 1:k
        for i in start[c]:(start[c + 1] - 1)
            v = members[i]
            for p in nzrange(A, v)
                d = membership[rows[p]]
                if d == c
                    internal[c] += vals[p]
                else
                    add!(acc, d, vals[p])
                end
            end
        end
        ids = active(acc)
        sort!(ids)
        for d in ids
            rowval[pos] = d
            nzval[pos] = acc.weight[d]
            pos += 1
        end
        reset!(acc)
        colptr[c + 1] = pos
    end
    resize!(rowval, pos - 1)
    resize!(nzval, pos - 1)
    internal ./= 2   # every internal edge was seen from both endpoints
    return SparseMatrixCSC(k, k, colptr, rowval, nzval), internal
end
