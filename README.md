# Fixing Collisions

* [Google Doc](https://docs.google.com/document/d/13MetSJoTTZ0ptT0SERbst1SgG-KbgK48hozhko6mJxc/edit?usp=sharing)
* [LaTeX Write-up](https://www.overleaf.com/6555952782nttqwfwgksjb)

## API

### Collision Detection

```c++
detect_edge_vertex_collisions(
    vertices, displacements, edges, method=(BRUTE_FORCE | HASH_MAP))
    // -> vector(vertex_idx, edge_idx, time_of_impact)
```

### Collision Pruning

```c++
prune_edge_vertex_collisions(vector(vertex_idx, edge_idx, time_of_impact))
    // -> vector(vertex_idx, edge_idx, time_of_impact)
```

### Space-Time Interference Volume Computation

```c++
compute_volume(
    vertices, displacements, edges,
    vector(vertex_idx, edge_idx, time_of_impact))
    // -> double(volume)
```

### Space-Time Interference Volume Gradient Computation

```c++
compute_volume_gradient(
    vertices, displacements, edges,
    vector(vertex_idx, edge_idx, time_of_impact))
    // -> vector([d volume / d displacement for displacement in displacements])
```
