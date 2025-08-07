# connected.jl

"""
    is_connected_community(g::AbstractGraph, nodes::Set{Int})

Check if the nodes form a connected subgraph.
"""
function is_connected_community(g::AbstractGraph, nodes::Set{Int})
    if isempty(nodes)
        return true
    end
    
    # Start BFS from arbitrary node
    start_node = first(nodes)
    visited = Set{Int}([start_node])
    queue = [start_node]
    
    while !isempty(queue)
        current = popfirst!(queue)
        for neighbor in neighbors(g, current)
            if neighbor in nodes && !(neighbor in visited)
                push!(visited, neighbor)
                push!(queue, neighbor)
            end
        end
    end
    
    return length(visited) == length(nodes)
end

"""
    find_connected_components(g::AbstractGraph, nodes::Set{Int})

Find all connected components within a set of nodes.
"""
function find_connected_components(g::AbstractGraph, nodes::Set{Int})
    components = Vector{Set{Int}}()
    unvisited = copy(nodes)
    
    while !isempty(unvisited)
        # Start new component
        start_node = first(unvisited)
        component = Set{Int}([start_node])
        delete!(unvisited, start_node)
        queue = [start_node]
        
        # BFS to find all connected nodes
        while !isempty(queue)
            current = popfirst!(queue)
            for neighbor in neighbors(g, current)
                if neighbor in unvisited
                    push!(component, neighbor)
                    delete!(unvisited, neighbor)
                    push!(queue, neighbor)
                end
            end
        end
        
        push!(components, component)
    end
    
    return components
end