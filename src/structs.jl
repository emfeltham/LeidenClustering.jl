# Leiden.jl

"""
    LeidenPartition

Represents a partition of nodes into communities.
"""
mutable struct LeidenPartition
    membership::Vector{Int}           # Node -> Community mapping
    community_nodes::Dict{Int, Set{Int}}  # Community -> Set of nodes
    modularity::Float64
    resolution::Float64              # Resolution parameter
end

"""
    LeidenState

Maintains the state of the Leiden algorithm during execution.
"""
mutable struct LeidenState
    graph::AbstractGraph
    weights::SparseMatrixCSC{Float64}
    node_degrees::Vector{Float64}    # Sum of edge weights for each node
    total_weight::Float64            # Sum of all edge weights
    partition::LeidenPartition
    community_weights_in::Dict{Int, Float64}  # Internal weight of each community
    community_weights_tot::Dict{Int, Float64} # Total weight of edges incident to community
end
