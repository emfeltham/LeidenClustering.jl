# API reference

## Clustering

```@docs
leiden_clustering
louvain_clustering
```

## Results

```@docs
Partition
ncommunities
```

## Quality

`modularity` is not exported (`Graphs.modularity` has the same name); call it as
`LeidenClustering.modularity`.

```@docs
LeidenClustering.modularity(::Graphs.AbstractGraph, ::AbstractVector{<:Integer})
```

## Module

```@docs
LeidenClustering
```
