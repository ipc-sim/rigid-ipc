# Bullet Comparison

Bullet provides two ways of modeling contacts between complex geometries. Similar to our method, the first way uses the triangle mesh geometry directly (`btGImpactMeshShape`). The second way uses convex shapes (or compound convex shapes). Bullet additionally, provides a direct way of generating approximate convex decompositions of concave shapes through Volumetric Hierarchical Approximate Decomposition (V-HACD).

## Triangle Mesh Collisions (`btGImpactMeshShape`)

We modified the ImportMJCFDemo example to automatically take an MJCF scene file as input, run the scene for 10 simulated seconds, and then terminate. We use `GImpact` for any input triangle meshes, and enforce only one substep for each time step (time step size fixed at `1e-2`, `1e-3`, or `1e-4` s for our experiments). To do this we implemented mesh volume computation and enabled automatic mass computation from density input (`m=ρV`), and we also enabled mesh scaling, `dt`, and `mu` setting from MJCF file which are not supported in the original demo.

See our repository on Github for more details: [https://github.com/liminchen/bullet3](https://github.com/liminchen/bullet3)

Scene MJCF files can be found in the `MJCF` directory (files with extension `.xml`).

## Convex Shapes and V-HACD

We additionally compare against the convex decompositions of our meshes. We do this directly using PyBullet Python bindings. We directly load our JSON fixtures and use V-HACD to decompose concave meshes.

### Prerequisites

PyBullet Python bindings can be installed using:
```
pip install pybullet
```

For videos, make sure `ffmpeg` is installed and in the path.

### Running Simulation

To run the comparison run `python benchmark.py`.
