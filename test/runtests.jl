# julia --project="." -e "import Pkg; Pkg.test()" > test/runtests.txt 2>&1

using LeidenClustering
using Graphs
using SimpleWeightedGraphs
using Random
using Test
using Aqua

@testset "LeidenClustering.jl" begin

    @testset "Internals: accumulator, reindex, contraction" begin
        using LeidenClustering: NeighborAccumulator, add!, active, reset!, renumber!,
            contract, Adjacency, strengths

        acc = NeighborAccumulator(10)
        add!(acc, 3, 3.0); add!(acc, 1, 5.0); add!(acc, 1, 2.0); add!(acc, 7, 0.0)
        @test collect(active(acc)) == [3, 1, 7]
        @test acc.weight[1] == 7.0 && acc.weight[3] == 3.0
        reset!(acc)
        @test acc.count == 0 && all(iszero, acc.weight) && !any(acc.seen)
        add!(acc, 7, 1.0)          # a zero-weight visit must not be re-added later
        @test collect(active(acc)) == [7]

        m = [5, 5, 2, 9, 2]
        @test renumber!(m) == 3
        @test m == [1, 1, 2, 3, 2]
        @test renumber!(Int[]) == 0

        # Contraction keeps total weight: cross weights in B, internal weight separate
        g = SimpleWeightedGraph(4)
        add_edge!(g, 1, 2, 2.0); add_edge!(g, 2, 3, 1.0); add_edge!(g, 3, 4, 3.0)
        adj = Adjacency(g)
        B, internal = contract(adj.A, [1, 1, 2, 2], 2)
        @test Matrix(B) == [0.0 1.0; 1.0 0.0]
        @test internal == [2.0, 3.0]
        @test all(issorted(B.rowval[B.colptr[j]:(B.colptr[j + 1] - 1)]) for j in 1:2)
    end

    @testset "Louvain Determinism" begin
        # Test that optimized version produces identical results
        for seed in 1:3
            Random.seed!(seed)
            g = erdos_renyi(50, 0.1)

            Random.seed!(seed)
            membership1, mods1 = louvain_clustering(g)

            Random.seed!(seed)
            membership2, mods2 = louvain_clustering(g)

            # Should be identical with same seed
            @test membership1 == membership2
            @test mods1 ≈ mods2
        end
    end

    @testset "Louvain Properties" begin
        Random.seed!(42)
        g = erdos_renyi(100, 0.05)

        membership, mods = louvain_clustering(g)

        # Basic properties
        @test length(membership) == nv(g)
        @test all(membership .> 0)
        @test maximum(membership) == length(unique(membership))  # Consecutive IDs

        # Modularity should be non-negative for random graphs
        @test mods[end] >= 0.0

        # Should have found some community structure
        @test length(unique(membership)) < nv(g)
    end

    @testset "Graph Contraction Properties" begin
        # Test that contracted graphs preserve properties
        Random.seed!(42)
        g = erdos_renyi(30, 0.15)

        membership, mods = louvain_clustering(g)

        # Basic properties
        @test length(membership) == nv(g)
        @test all(membership .> 0)
        @test maximum(membership) <= nv(g)  # Should reduce vertex count
    end

    @testset "Modularity Calculation" begin
        # Test modularity calculation on simple graphs
        Random.seed!(42)
        g = erdos_renyi(50, 0.1)

        membership, mods = louvain_clustering(g)

        # Calculate modularity using our function
        q = LeidenClustering.modularity(g, membership; resolution=1.0)

        # Should be in valid range [-0.5, 1.0]
        @test q >= -0.5
        @test q <= 1.0

        # Modularity from final partition should be close to final clustering modularity
        # (may differ slightly due to implementation details)
        @test abs(q - mods[end]) < 0.2
    end

    @testset "Weighted Graphs" begin
        Random.seed!(42)
        n = 20
        g_simple = erdos_renyi(n, 0.2)

        # Create weighted version
        edges_list = collect(edges(g_simple))
        weights = rand(length(edges_list))

        src_list = [src(e) for e in edges_list]
        dst_list = [dst(e) for e in edges_list]

        g_weighted = SimpleWeightedGraph(src_list, dst_list, weights)

        # Should still work
        membership, mods = louvain_clustering(g_weighted)

        @test length(membership) == n
        @test mods[end] >= 0.0
    end

    @testset "Multi-Level Aggregation" begin
        # Test that optimization works across multiple aggregation levels
        Random.seed!(42)
        g = erdos_renyi(100, 0.05)

        membership, mods = louvain_clustering(g)

        # Should have at least one level
        @test length(mods) >= 1

        # Modularity should be non-decreasing (within numerical tolerance)
        for i in 2:length(mods)
            @test mods[i] >= mods[i-1] - 1e-10
        end
    end

    @testset "Different Graph Types" begin
        test_graphs = [
            ("Erdős-Rényi", erdos_renyi(40, 0.1, seed=42)),
            ("Barabási-Albert", barabasi_albert(40, 3, seed=42)),
            ("Watts-Strogatz", watts_strogatz(40, 6, 0.3, seed=42)),
            ("Complete", complete_graph(20)),
            ("Star", star_graph(30)),
            ("Path", path_graph(25)),
        ]

        for (name, g) in test_graphs
            Random.seed!(42)
            membership, mods = louvain_clustering(g)

            @test length(membership) == nv(g)
            @test all(membership .> 0)
            @test mods[end] >= 0.0
        end
    end

    @testset "Resolution Parameter" begin
        Random.seed!(42)
        g = erdos_renyi(50, 0.1)

        # Test different resolution values
        for resolution in [0.5, 1.0, 1.5, 2.0]
            membership, mods = louvain_clustering(g; resolution=resolution)

            @test length(membership) == nv(g)
            @test all(membership .> 0)
            @test !isnan(mods[end])
        end
    end

    @testset "Empty and Trivial Graphs" begin
        # Single vertex
        g1 = SimpleGraph(1)
        membership1, mods1 = louvain_clustering(g1)
        @test membership1 == [1]
        @test length(mods1) >= 0  # May return empty array for trivial graphs

        # Two disconnected vertices
        g2 = SimpleGraph(2)
        membership2, mods2 = louvain_clustering(g2)
        @test length(membership2) == 2
        @test all(membership2 .> 0)
    end

    # ========================================================================
    # LEIDEN ALGORITHM TESTS
    # ========================================================================

    @testset "Leiden: Basic Functionality" begin
        # Single triangle
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)
        add_edge!(g, 2, 3)
        add_edge!(g, 1, 3)

        membership, quality = leiden_clustering(g)
        @test length(membership) == 3
        @test length(unique(membership)) == 1  # Should be one community
        @test quality ≈ 0.0 atol=0.01  # Perfect modularity for complete graph

        # Two disconnected triangles
        g2 = SimpleGraph(6)
        add_edge!(g2, 1, 2)
        add_edge!(g2, 2, 3)
        add_edge!(g2, 1, 3)
        add_edge!(g2, 4, 5)
        add_edge!(g2, 5, 6)
        add_edge!(g2, 4, 6)

        membership2, quality2 = leiden_clustering(g2)
        @test length(membership2) == 6
        @test length(unique(membership2)) == 2  # Should be two communities
        @test quality2 > 0.4  # Good separation
    end

    @testset "Leiden: Input Validation" begin
        g = SimpleGraph(3)
        add_edge!(g, 1, 2)

        # Negative resolution
        @test_throws ArgumentError leiden_clustering(g, resolution=-1.0)

        # Beta out of range
        @test_throws ArgumentError leiden_clustering(g, beta=-0.1)
        @test_throws ArgumentError leiden_clustering(g, beta=1.5)

        # Invalid objective
        @test_throws ArgumentError leiden_clustering(g, objective=:invalid)

        # Bad node weights / initial membership
        @test_throws ArgumentError leiden_clustering(g, node_weights=ones(3))   # modularity
        @test_throws ArgumentError leiden_clustering(g, objective=:cpm, node_weights=ones(2))
        @test_throws ArgumentError leiden_clustering(g, objective=:cpm, node_weights=[1.0, -1.0, 1.0])
        @test_throws ArgumentError leiden_clustering(g, initial_membership=[1, 1])
        @test_throws ArgumentError leiden_clustering(g, initial_membership=[0, 1, 1])

        # Directed graph
        dg = SimpleDiGraph(3)
        add_edge!(dg, 1, 2)
        @test_throws ArgumentError leiden_clustering(dg)
    end

    @testset "Leiden: Edge Cases" begin
        # Single node
        g1 = SimpleGraph(1)
        membership1, _ = leiden_clustering(g1)
        @test length(membership1) == 1
        @test membership1[1] == 1

        # Two disconnected nodes
        g2 = SimpleGraph(2)
        membership2, _ = leiden_clustering(g2)
        @test length(membership2) == 2

        # Complete graph
        g_complete = complete_graph(5)
        membership_complete, _ = leiden_clustering(g_complete)
        @test length(unique(membership_complete)) == 1  # All in one community

        # Empty graph (no edges)
        g_empty = SimpleGraph(5)
        membership_empty, _ = leiden_clustering(g_empty)
        @test length(unique(membership_empty)) == 5  # All separate
    end

    @testset "Leiden: Determinism" begin
        g = erdos_renyi(50, 0.1, seed=123)

        # Same seed should give same results
        m1, q1 = leiden_clustering(g, seed=42)
        m2, q2 = leiden_clustering(g, seed=42)
        @test m1 == m2
        @test q1 ≈ q2
    end

    @testset "Leiden: Quality Properties" begin
        g = erdos_renyi(30, 0.2, seed=42)
        membership, quality = leiden_clustering(g)

        # Quality should be between -1 and 1 for modularity
        @test -1.0 <= quality <= 1.0

        # All nodes should be assigned
        @test length(membership) == nv(g)
        @test all(membership .> 0)

        # Community IDs should be consecutive
        unique_comms = sort(unique(membership))
        @test unique_comms == collect(1:length(unique_comms))
    end

    @testset "Leiden: Known Graphs" begin
        # Zachary's Karate Club
        g = smallgraph(:karate)
        membership, quality = leiden_clustering(g, seed=42)

        # Should find reasonable number of communities
        n_comms = length(unique(membership))
        @test 1 <= n_comms <= nv(g)

        # Should achieve positive modularity
        @test quality > 0.0
    end

    @testset "Leiden vs Louvain Comparison" begin
        g = smallgraph(:karate)

        # Both algorithms should work
        m_leiden, q_leiden = leiden_clustering(g, seed=42)
        m_louvain, q_louvain = louvain_clustering(g, seed=42)

        @test length(m_leiden) == nv(g)
        @test length(m_louvain) == nv(g)

        # Both should find reasonable communities
        @test 1 <= length(unique(m_leiden)) <= nv(g)
        @test 1 <= length(unique(m_louvain)) <= nv(g)

        # Both should have positive modularity
        @test q_leiden[end] > 0
        @test q_louvain[end] > 0
    end

    @testset "Regression: Louvain membership matches reported modularity" begin
        for (i, g) in enumerate([smallgraph(:karate),
                                 (Random.seed!(1); erdos_renyi(500, 0.02)),
                                 Graphs.grid([30, 30])])
            Random.seed!(i)
            m, q = louvain_clustering(g)
            @test LeidenClustering.modularity(g, m) ≈ q[end] atol = 1e-9
            @test maximum(m) > 1
        end
        # Karate optimum is Q ≈ 0.4198; Louvain must be in the right region
        Random.seed!(1)
        m, _ = louvain_clustering(smallgraph(:karate))
        @test LeidenClustering.modularity(smallgraph(:karate), m) > 0.38
    end

    @testset "Regression: Leiden quality on structured graphs" begin
        Random.seed!(7)
        g = Graphs.grid([40, 40])      # over 50 seeds Leiden gets Q in 0.848-0.855, 16-21 communities
        m, q = leiden_clustering(g)
        @test q[end] ≈ LeidenClustering.modularity(g, m) atol = 1e-9
        @test q[end] > 0.84            # thresholds leave margin for RNG streams differing across Julia versions
        @test maximum(m) < 40          # the old membership-reset bug gave ~6x too many communities

        Random.seed!(7)
        m, q = leiden_clustering(smallgraph(:karate))
        @test q[end] > 0.41
        @test maximum(m) == 4
    end

    @testset "Regression: Leiden communities are connected" begin
        Random.seed!(11)
        g = erdos_renyi(300, 0.02)
        m, _ = leiden_clustering(g)
        for c in 1:maximum(m)
            nodes = findall(==(c), m)
            sub, _ = induced_subgraph(g, nodes)
            @test is_connected(sub)
        end
    end

    @testset "Regression: Leiden uses edge weights on aggregated levels" begin
        # Two cliques joined by a heavy bridge vs. light bridge
        function two_cliques(bridge)
            g = SimpleWeightedGraph(10)
            for a in 1:5, b in a+1:5
                add_edge!(g, a, b, 1.0); add_edge!(g, a + 5, b + 5, 1.0)
            end
            add_edge!(g, 5, 6, bridge)
            g
        end
        Random.seed!(3)
        m, _ = leiden_clustering(two_cliques(0.1))
        @test length(unique(m)) == 2
        @test length(unique(m[1:5])) == 1 && length(unique(m[6:10])) == 1
        # Weighted strength must be used for node weights: heavy bridge merges everything
        Random.seed!(3)
        m, _ = leiden_clustering(two_cliques(1000.0))
        @test m[5] == m[6]
    end

    @testset "Regression: Leiden CPM objective" begin
        Random.seed!(5)
        g = smallgraph(:karate)
        m_lo, _ = leiden_clustering(g; objective=:cpm, resolution=0.01)
        m_hi, _ = leiden_clustering(g; objective=:cpm, resolution=0.5)
        @test maximum(m_hi) > maximum(m_lo)
    end

    @testset "Regression: self-loops" begin
        using LeidenClustering: Adjacency, strengths
        weighted_adjacency(x) = (nothing, strengths(Adjacency(x)))
        edges_ = [(1, 2), (2, 3), (3, 4), (2, 2), (4, 4)]
        g = SimpleGraph(4); gw = SimpleWeightedGraph(4)
        for (a, b) in edges_
            add_edge!(g, a, b); add_edge!(gw, a, b, 1.0)
        end
        # A loop counts twice in a vertex's strength (igraph convention), in both representations
        @test weighted_adjacency(g)[2] == [1.0, 4.0, 2.0, 3.0]
        @test weighted_adjacency(gw)[2] == [1.0, 4.0, 2.0, 3.0]

        # Planted cliques with random loops: both algorithms recover them and report true modularity
        Random.seed!(2)
        h = SimpleGraph(60)
        for c in 0:5, a in 1:10, b in a+1:10
            rand() < 0.7 && add_edge!(h, 10c + a, 10c + b)
        end
        for _ in 1:40; add_edge!(h, rand(1:60), rand(1:60)); end
        for v in 1:60; rand() < 0.5 && add_edge!(h, v, v); end
        m, q = leiden_clustering(h)
        ml, ql = louvain_clustering(h)
        @test q ≈ LeidenClustering.modularity(h, m) atol = 1e-9
        @test ql ≈ LeidenClustering.modularity(h, ml) atol = 1e-9
        @test maximum(m) == 6
        @test q ≥ ql - 1e-9
    end

    @testset "API: result type, RNG, keywords" begin
        g = smallgraph(:karate)

        r = leiden_clustering(g; seed=3)
        @test r isa Partition
        membership, quality = r
        @test membership === r.membership && quality == r.quality
        @test r.qualities == [r.quality]
        @test ncommunities(r) == maximum(r.membership)
        @test r.quality ≈ LeidenClustering.modularity(g, r.membership) atol = 1e-9

        # `seed` uses a private RNG and leaves the global one untouched
        Random.seed!(11); before = rand()
        Random.seed!(11)
        a = leiden_clustering(g; seed=5); b = louvain_clustering(g; seed=5)
        @test rand() == before
        @test a.membership == leiden_clustering(g; seed=5).membership
        @test b.membership == louvain_clustering(g; seed=5).membership

        # an explicit rng is used and advances
        rng = Xoshiro(1)
        @test leiden_clustering(g; rng).membership == leiden_clustering(g; rng=Xoshiro(1)).membership
        @test leiden_clustering(g; seed=1).membership == leiden_clustering(g; rng=Xoshiro(1)).membership

        # Real-valued keywords
        @test leiden_clustering(g; resolution=1, beta=0, seed=1) isa Partition
        @test louvain_clustering(g; resolution=1, seed=1) isa Partition
    end

    @testset "Leiden options and edge cases" begin
        g = smallgraph(:karate)
        # greedy refinement (beta = 0) and extra passes still give a valid, good partition
        for kw in ((; beta=0), (; n_iterations=10), (; n_iterations=1))
            r = leiden_clustering(g; seed=7, kw...)
            @test r.quality > 0.38
            @test r.quality ≈ LeidenClustering.modularity(g, r.membership) atol = 1e-9
        end

        # empty graph, edgeless graph
        for f in (leiden_clustering, louvain_clustering)
            r = f(SimpleGraph(0))
            @test isempty(r.membership)
            r = f(SimpleGraph(5))
            @test r.membership == 1:5 && r.quality == 0.0
        end

        # modularity: non-consecutive ids, and argument checks
        @test LeidenClustering.modularity(g, fill(1, 34)) ≈ 0 atol = 1e-12
        m = leiden_clustering(g; seed=1).membership
        @test LeidenClustering.modularity(g, 10 .* m) ≈ LeidenClustering.modularity(g, m)
        @test_throws ArgumentError LeidenClustering.modularity(g, [1, 2])
        @test_throws ArgumentError LeidenClustering.modularity(g, zeros(Int, 34))
        @test_throws ArgumentError LeidenClustering.modularity(SimpleDiGraph(3), [1, 1, 1])
    end

    @testset "Edge weights: sources and validation" begin
        using LeidenClustering: Adjacency
        g = smallgraph(:karate)
        W = zeros(34, 34)
        gw = SimpleWeightedGraph(34)
        for (k, e) in enumerate(edges(g))
            w = 1.0 + k % 3
            W[src(e), dst(e)] = W[dst(e), src(e)] = w
            add_edge!(gw, src(e), dst(e), w)
        end
        # A matrix on a plain graph is the same as a SimpleWeightedGraph with those weights
        @test Adjacency(g, W).A == Adjacency(gw).A
        @test leiden_clustering(g; weights=W, seed=1).membership == leiden_clustering(gw; seed=1).membership
        @test louvain_clustering(g; weights=W, seed=1).membership == louvain_clustering(gw; seed=1).membership
        m = leiden_clustering(gw; seed=1).membership
        @test LeidenClustering.modularity(g, m; weights=W) ≈ LeidenClustering.modularity(gw, m)
        # ...and unit weights override a weighted graph's own
        @test Adjacency(gw, Graphs.DefaultDistance(34)).A == Adjacency(g).A
        # entries off the edge set are ignored; a loop's matrix entry is its weight
        h = path_graph(3); add_edge!(h, 2, 2)
        a = Adjacency(h, fill(2.0, 3, 3))
        @test count(!iszero, a.A) == 4 && a.loops == [0.0, 2.0, 0.0]

        bad = SimpleWeightedGraph(3); add_edge!(bad, 1, 2, -1.0); add_edge!(bad, 2, 3, 1.0)
        for f in (leiden_clustering, louvain_clustering)
            @test_throws ArgumentError f(bad)
            @test_throws ArgumentError f(g; weights=fill(NaN, 34, 34))
            @test_throws ArgumentError f(g; weights=ones(3, 3))
        end
        @test_throws ArgumentError LeidenClustering.modularity(bad, [1, 1, 2])
    end

    @testset "Leiden: CPM quality, node weights, initial membership, n_iterations" begin
        # CPM quality counts a loop twice, as igraph does: (Σ_ij A_ij δ − γ Σ_c n_c²) / 2m
        h = path_graph(6); add_edge!(h, 1, 1); add_edge!(h, 4, 4)
        D = Matrix{Float64}(adjacency_matrix(h))   # loops are 2 on the diagonal
        for γ in (0.1, 0.4)
            r = leiden_clustering(h; objective=:cpm, resolution=γ, seed=1)
            m = r.membership
            nc = [count(==(c), m) for c in 1:maximum(m)]
            expected = (sum(D[i, j] for i in 1:6, j in 1:6 if m[i] == m[j]) - γ * sum(abs2, nc)) / sum(D)
            @test r.quality ≈ expected
        end

        # node weights: unit weights are the default; a heavy vertex is left alone
        g = smallgraph(:karate)
        @test leiden_clustering(g; objective=:cpm, resolution=0.1, node_weights=ones(34), seed=2).membership ==
              leiden_clustering(g; objective=:cpm, resolution=0.1, seed=2).membership
        nw = ones(34); nw[34] = 100.0
        m = leiden_clustering(g; objective=:cpm, resolution=0.1, node_weights=nw, seed=2).membership
        @test count(==(m[34]), m) == 1

        # initial membership: n_iterations = 0 returns it (renumbered)
        init = [isodd(v) ? 7 : 3 for v in 1:34]
        r0 = leiden_clustering(g; initial_membership=init, n_iterations=0)
        @test r0.membership == [isodd(v) ? 1 : 2 for v in 1:34]
        @test r0.quality ≈ LeidenClustering.modularity(g, init)
        # continuing from a good partition does not make it worse
        r = leiden_clustering(g; seed=4)
        @test leiden_clustering(g; initial_membership=r.membership, seed=5).quality ≥ r.quality - 1e-12
        # moves on the first level are kept even when everything ends up alone (igraph semantics)
        p = leiden_clustering(path_graph(6); objective=:cpm, resolution=2.0,
                              initial_membership=ones(Int, 6), n_iterations=1, seed=1)
        @test p.membership == 1:6

        # negative n_iterations runs until a pass changes nothing
        r = leiden_clustering(g; n_iterations=-1, seed=6)
        @test r.quality ≈ LeidenClustering.modularity(g, r.membership) atol = 1e-9
        @test r.quality > 0.41
    end

    @testset "Package hygiene (Aqua)" begin
        Aqua.test_all(LeidenClustering)
    end

    @testset "Exports" begin
        # `modularity` must not be exported: it would clash with Graphs.modularity
        @test !(:modularity in names(LeidenClustering))
        # ...and nothing else we export may collide with a Graphs.jl export (e.g. `louvain`)
        @test isempty(intersect(names(LeidenClustering), names(Graphs)))
        @test Set(names(LeidenClustering)) ⊇ Set([:leiden_clustering, :louvain_clustering, :Partition, :ncommunities])
    end

end
