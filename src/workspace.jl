# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Eric Martin Feltham

# Small reusable building blocks for the local-moving loops.

"""
    NeighborAccumulator(n)

Sparse accumulator of weight per community id in `1:n`. `add!` sums weight into a
community, remembering which ones were used, so `reset!` costs O(#used) instead
of O(n). Iterate the used communities with [`active`](@ref).
"""
mutable struct NeighborAccumulator
    weight::Vector{Float64}
    seen::Vector{Bool}
    ids::Vector{Int}
    count::Int
end

NeighborAccumulator(n::Integer) =
    NeighborAccumulator(zeros(Float64, n), fill(false, n), Vector{Int}(undef, n), 0)

@inline function add!(acc::NeighborAccumulator, c::Int, w::Float64)
    @inbounds begin
        if !acc.seen[c]
            acc.seen[c] = true
            acc.count += 1
            acc.ids[acc.count] = c
        end
        acc.weight[c] += w
    end
    return acc
end

"Communities that received weight since the last `reset!`, in order of first appearance."
@inline active(acc::NeighborAccumulator) = @inbounds view(acc.ids, 1:acc.count)

function reset!(acc::NeighborAccumulator)
    @inbounds for i in 1:acc.count
        c = acc.ids[i]
        acc.weight[c] = 0.0
        acc.seen[c] = false
    end
    acc.count = 0
    return acc
end

# `rng` wins unless a `seed` is given, which creates a private generator: the global
# RNG is never reseeded.
resolve_rng(rng::AbstractRNG, seed::Nothing) = rng
resolve_rng(::AbstractRNG, seed::Integer) = Xoshiro(seed)
