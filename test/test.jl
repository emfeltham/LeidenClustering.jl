# testing.jl

import CSV
using DataFrames
using Graphs
using LeidenClustering

# Run test
g = make_simplegraph()
test_leiden(g)

karate = karateclub_graph();

test_leiden(karate)

test_resolution_sweep(karate, resolutions=[0.3, 0.5, 0.7, 1.0, 1.5])

test_resolution_sweep(karate, resolutions=[1.0])

#=
r$> cat("Number of communities:", length(sizes), "\n")
Number of communities: 4 

r$> print(sizes)
Community sizes
 1  2  3  4 
11  5  6 12 

r$> cat("Modularity:", mod, "\n")
Modularity: 0.4197896 
=#
