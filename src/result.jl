# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2025 Eric Martin Feltham

"""
    Partition

Result of [`leiden_clustering`](@ref) or [`louvain_clustering`](@ref).

# Fields
- `membership::Vector{Int}`: community id `1:k` of every vertex.
- `quality::Float64`: final quality of `membership` (modularity, or CPM for `objective=:cpm`).
- `qualities::Vector{Float64}`: quality after each level for Louvain; `[quality]` for Leiden.

A `Partition` destructures as `membership, quality = result`.
"""
struct Partition
    membership::Vector{Int}
    quality::Float64
    qualities::Vector{Float64}
end

Base.iterate(c::Partition, state::Int=1) =
    state == 1 ? (c.membership, 2) : state == 2 ? (c.quality, 3) : nothing

"Number of communities in `c`."
ncommunities(c::Partition) = maximum(c.membership; init=0)

function Base.show(io::IO, c::Partition)
    print(io, "Partition(", length(c.membership), " vertices, ", ncommunities(c),
          " communities, quality = ", round(c.quality; digits=4), ")")
end
