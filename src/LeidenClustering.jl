module LeidenClustering

using Graphs, SimpleWeightedGraphs
using SparseArrays
using Random

# example
using DataFrames: DataFrame
import CSV

include("structs.jl")

include("connected.jl")
include("local_moving.jl")

include("mechanics.jl")
include("refine.jl")
include("leiden.jl")
export leiden
include("show.jl")
include("testing.jl")
export
    make_simplegraph,
    test_leiden,
    test_resolution_sweep

end # module
