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
function test_leiden(g::AbstractGraph = nothing)
    if g === nothing
        # Create test graph: two triangles with bridge
        g = SimpleGraph(6)
        add_edge!(g, 1, 2)
        add_edge!(g, 1, 3)
        add_edge!(g, 2, 3)
        add_edge!(g, 4, 5)
        add_edge!(g, 4, 6)
        add_edge!(g, 5, 6)
        add_edge!(g, 3, 4)
    end
    
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
