# testing.jl

function make_simplegraph()
    # Create a simple test graph (karate club or similar)
    g = SimpleGraph(6)
    add_edge!(g, 1, 2)
    add_edge!(g, 1, 3)
    add_edge!(g, 2, 3)
    add_edge!(g, 4, 5)
    add_edge!(g, 4, 6)
    add_edge!(g, 5, 6)
    add_edge!(g, 3, 4)  # Bridge between communities
    return g
end

"""
    make_weighted_simplegraph()

Create a simple weighted test graph with community structure.
Requires SimpleWeightedGraphs to be loaded.
"""
function make_weighted_simplegraph()
    if !isdefined(Main, :SimpleWeightedGraph)
        error("SimpleWeightedGraphs package must be loaded to create weighted graphs")
    end
    
    g = Main.SimpleWeightedGraph(6)
    
    # Community 1: nodes 1,2,3 with strong internal connections
    add_edge!(g, 1, 2, 5.0)
    add_edge!(g, 1, 3, 5.0)
    add_edge!(g, 2, 3, 5.0)
    
    # Community 2: nodes 4,5,6 with strong internal connections
    add_edge!(g, 4, 5, 5.0)
    add_edge!(g, 4, 6, 5.0)
    add_edge!(g, 5, 6, 5.0)
    
    # Weak bridge between communities
    add_edge!(g, 3, 4, 0.1)
    
    return g
end

function karateclub_graph()
    kn = CSV.read("data/karate.csv", DataFrame);
    n = sort(unique(vcat(kn.src, kn.dst))) |> length;

    karate = SimpleGraph(n);
    for (a, b) in zip(kn.src, kn.dst)
        add_edge!(karate, a, b)
    end
    return karate
end

# Sanity check function
function sanity_check(state::LeidenState)
    """
    Check for common issues in the state.
    """
    println("\nSanity Check:")
    
    # Check 1: All nodes assigned
    n = nv(state.graph)
    assigned_nodes = Set{Int}()
    for (comm, nodes) in state.partition.community_nodes
        union!(assigned_nodes, nodes)
    end
    println("  All nodes assigned: $(length(assigned_nodes) == n) ($(length(assigned_nodes))/$n)")
    
    # Check 2: Membership consistency
    membership_consistent = true
    for (comm, nodes) in state.partition.community_nodes
        for node in nodes
            if state.partition.membership[node] != comm
                membership_consistent = false
                println("    ERROR: Node $node in community $comm but membership says $(state.partition.membership[node])")
            end
        end
    end
    println("  Membership consistent: $membership_consistent")
    
    # Check 3: No -1 memberships
    has_neg_one = any(x -> x == -1, state.partition.membership)
    println("  No -1 memberships: $(!has_neg_one)")
    
    # Check 4: Weight consistency
    total_degree = sum(state.node_degrees)
    comm_total = sum(values(state.community_weights_tot))
    println("  Total degree match: $(abs(total_degree - comm_total) < 1e-10) (diff: $(abs(total_degree - comm_total)))")
    
    # Check 5: Internal weights reasonable
    total_internal = sum(values(state.community_weights_in))
    println("  Total internal weight: $total_internal (should be ≤ $(state.total_weight))")
    
    return membership_consistent && !has_neg_one
end

# Testing functions
function test_leiden(g)

    println("\nTesting Full Leiden Algorithm")
    println("="^50)
    println("Graph: $(nv(g)) nodes, $(ne(g)) edges")
    
    # Test without refinement (pure Louvain)
    println("\n1. Pure Louvain (no refinement):")
    state_louvain = leiden(g, resolution=1.0, use_refinement=false, seed=42)
    println("  Communities: $(length(state_louvain.partition.community_nodes))")
    println("  Modularity: $(round(state_louvain.partition.modularity, digits=4))")
    
    # Test with refinement (Leiden)
    println("\n2. Leiden (with refinement):")
    state_leiden = leiden(g, resolution=1.0, use_refinement=true, seed=42)
    println("  Communities: $(length(state_leiden.partition.community_nodes))")
    println("  Modularity: $(round(state_leiden.partition.modularity, digits=4))")
    
    return state_louvain, state_leiden
end

function test_resolution_sweep(g::AbstractGraph; resolutions=[0.5, 1.0, 1.5], seed=42)
    println("\nResolution Parameter Sweep (with full Leiden)")
    println("="^50)
    
    for res in resolutions
        state = leiden(g, resolution=res, seed=seed)
        
        # Count community sizes
        sizes = sort([length(nodes) for (_, nodes) in state.partition.community_nodes], rev=true)
        
        println("\nResolution: $res")
        println("  Modularity: $(round(state.partition.modularity, digits=4))")
        println("  Communities: $(length(state.partition.community_nodes))")
        println("  Sizes: $sizes")
    end
end

# Diagnostic function to verify modularity calculation
function verify_modularity(state::LeidenState)
    """
    Manually calculate modularity to verify it's correct.
    """
    n = nv(state.graph)
    m = state.total_weight
    Q = 0.0
    
    for i in 1:n
        for j in 1:n
            if state.partition.membership[i] == state.partition.membership[j]
                A_ij = state.weights[i, j]
                k_i = state.node_degrees[i]
                k_j = state.node_degrees[j]
                Q += A_ij - state.partition.resolution * (k_i * k_j) / (2.0 * m)
            end
        end
    end
    
    Q = Q / (2.0 * m)
    
    println("\nModularity verification:")
    println("  Computed: $(round(state.partition.modularity, digits=6))")
    println("  Manual:   $(round(Q, digits=6))")
    println("  Match: $(abs(state.partition.modularity - Q) < 1e-6)")
    
    return Q
end


"""
    test_weighted_leiden()

Test Leiden algorithm on weighted graphs.
"""
function test_weighted_leiden()
    if !isdefined(Main, :SimpleWeightedGraph)
        println("SimpleWeightedGraphs not loaded, skipping weighted tests")
        return nothing
    end
    
    println("\nTesting Weighted Leiden Algorithm")
    println("="^50)
    
    # Create weighted graph
    g_weighted = make_weighted_simplegraph()
    println("Weighted graph: $(nv(g_weighted)) nodes, $(ne(g_weighted)) edges")
    
    # Test with weighted graph
    println("\n1. Weighted graph clustering:")
    state_weighted = leiden(g_weighted, resolution=1.0, seed=42)
    println("  Communities: $(length(state_weighted.partition.community_nodes))")
    println("  Modularity: $(round(state_weighted.partition.modularity, digits=4))")
    
    # Create unweighted version of the same topology
    g_unweighted = SimpleGraph(6)
    for e in edges(g_weighted)
        add_edge!(g_unweighted, src(e), dst(e))
    end
    
    println("\n2. Unweighted graph clustering (same topology):")
    state_unweighted = leiden(g_unweighted, resolution=1.0, seed=42)
    println("  Communities: $(length(state_unweighted.partition.community_nodes))")
    println("  Modularity: $(round(state_unweighted.partition.modularity, digits=4))")
    
    # The weighted version should separate communities better due to weak bridge
    println("\n3. Comparison:")
    println("  Weighted communities should separate at weak bridge (0.1 weight)")
    println("  Unweighted treats all edges equally (weight 1.0)")
    
    # Show the actual communities
    println("\n4. Community assignments:")
    println("  Weighted: $(state_weighted.partition.membership)")
    println("  Unweighted: $(state_unweighted.partition.membership)")
    
    return state_weighted, state_unweighted
end

"""
    test_weight_validation()

Test that negative weights are properly rejected.
"""
function test_weight_validation()
    println("\nTesting Weight Validation")
    println("="^50)
    
    g = SimpleGraph(3)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    
    # Create weight matrix with negative weight
    weights = sparse([1, 2, 2, 3], [2, 1, 3, 2], [1.0, 1.0, -1.0, -1.0], 3, 3)
    
    try
        state = initialize_state(g, weights=weights)
        println("ERROR: Should have rejected negative weights!")
    catch e
        if isa(e, ArgumentError)
            println("✓ Correctly rejected negative weights")
            println("  Error message: $(e.msg)")
        else
            println("ERROR: Wrong exception type: $(typeof(e))")
        end
    end
end