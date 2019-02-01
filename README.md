# Fixing Collisions

* [Google Doc](https://docs.google.com/document/d/13MetSJoTTZ0ptT0SERbst1SgG-KbgK48hozhko6mJxc/edit?usp=sharing)
* [LaTeX Write-up](https://www.overleaf.com/6555952782nttqwfwgksjb)

## API

### Collision Detection

```c++
detect_edge_vertex_collisions(vertices_at_t0, vertices_at_t1, edges); // -> vector(vertex_idx, edge_idx, time_of_impact)
detect_edge_vertex_collisions_brute_force(vertices_at_t0, vertices_at_t1, edges); // -> vector(vertex_idx, edge_idx, time_of_impact)
detect_edge_vertex_collisions_hash_map(vertices_at_t0, vertices_at_t1, edges); // -> vector(vertex_idx, edge_idx, time_of_impact)
```

### Collision Pruning 

```c++
prune_edge_vertex_collisions(vector(vertex_idx, edge_idx, time_of_impact)); // -> vector(vertex_idx, edge_idx, time_of_impact)
```

### Space-Time Interference Volume Computation

```c++
compute_volume(vertices_at_t0, vertices_at_t1, edges, vector(vertex_idx, edge_idx, time_of_impact)); // -> double(volume)
```

### Space-Time Interference Volume Gradient Computation

```c++
compute_volume_gradient(vertices_at_t0, vertices_at_t1, edges, vector(vertex_idx, edge_idx, time_of_impact)); // -> vector(d volume / d displacement_0_x, d volume / d displacement_0_y, d volume / d displacement_1_x, d volume / d displacement_1_y, d volume / d displacement_2_x, d volume / d displacement_2_y)
```
