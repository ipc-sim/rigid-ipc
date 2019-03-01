#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
from pathlib import Path

import numpy as np
import meshio
import json

from scenes.edge_topology import edge_topology

base_path = Path(os.path.dirname(os.path.realpath(__file__)))
fixtures_path = base_path.joinpath("..","..","fixtures")

def main(args=None):
    parser = argparse.ArgumentParser(args)
    parser.add_argument("output", type=str, help="path to the output folder")
    parser.add_argument("--num_links","-n", type=int, default=2)
    parser.add_argument("--dt", "-d", type=float, default=1e-2)
    parser.add_argument("--time", "-t", type=float, nargs=2, default=[45,46])
    args = parser.parse_args()

    num_links = args.num_links
    output = args.output
    t0, t1 = args.time

    dt = args.dt
    gravity = np.array([0.0, -9.81], dtype=np.float64)
    delta = gravity / 2.0 * (dt ** 2)

    base_link = meshio.read(str(fixtures_path.joinpath("meshes", "chain_link.vtk")))

    V = base_link.points.copy()
    F = base_link.cells['triangle'].copy()
    EV, FE, EF = edge_topology(F)
    boundary_edges = np.logical_or(EF[:,0] == -1 , EF[:,1] == -1)
    EVb = EV[boundary_edges, :]
    V = V[:,0:2]

    vertices = np.empty((0, 2), dtype=np.float64)
    displ = np.empty((0, 2), dtype=np.float64)
    edges = np.empty((0, 2), dtype=np.int32)

    def add_link(x, y, dynamic):
        nonlocal vertices
        nonlocal edges
        nonlocal displ

        if dynamic:
            delta0 = gravity / 2.0 * ((t0 * dt) ** 2)
            delta1 = gravity / 2.0 * ((t1 * dt) ** 2) - delta0
        else:
            delta0 = [0.0,0.0]
            delta1 = [0.0,0.0]

        new_vertices = V + [x,y]
        new_vertices += delta0
        edges = np.append(edges, EVb + vertices.shape[0], axis=0)
        vertices = np.append(vertices, new_vertices, axis=0)
        displ = np.append(displ, np.ones_like(V) * delta1, axis=0)

    y = 40.0
    for i in range(0, num_links):
        add_link(0.0, y, (True if i > 0 else False))
        y -= 9.0

    # scale to fit in our window 10 x 10 window
    all_nodes = np.append(vertices, vertices + displ, axis=0)
    bbox_min = np.min(all_nodes, axis=0)
    bbox_max = np.max(all_nodes, axis=0)
    center = np.average(all_nodes, axis=0)
    scale = np.linalg.norm(bbox_max - bbox_min) # diagonal

    vertices = (vertices - center) * 10 / scale
    displ *= 10 / scale

    data = dict(
        vertices=vertices.tolist(),
        edges =edges.tolist(),
        displacements=displ.tolist()
    )
    with open(output, 'w') as outfile:
        json.dump(data, outfile, indent=4)



if __name__ == "__main__":
    main()
