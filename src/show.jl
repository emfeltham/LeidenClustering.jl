# show.jl

# Function to show detailed community info
function show_communities(state::LeidenState)
    println("\nCommunity Details:")
    for (comm, nodes) in state.partition.community_nodes
        nodes_list = sort(collect(nodes))
        internal = state.community_weights_in[comm]
        total = state.community_weights_tot[comm]
        println("  Community $comm: $(length(nodes)) nodes")
        println("    Nodes: $(nodes_list[1:min(10, length(nodes_list))])$(length(nodes_list) > 10 ? "..." : "")")
        println("    Internal weight: $(round(internal, digits=2))")
        println("    Total weight: $(round(total, digits=2))")
    end
end