#!/usr/bin/env python

import argparse
import os
from pathlib import Path

import numpy as np
import meshio
import json

from fixture_utils import *


def link(a, b, c, d, e):
    #           e
    #         _____
    #       b|     |c
    # a|_____|  |__'__
    #  |  d  |  |  ,
    #        |_____|
    #            f

    V = np.array([[0.0, -a], [0.0, 0.0], [0.0, a], [d, -b], [d, 0.0], [d, b],
                  [d + e, -b], [d + e, -b + c], [d + e, +b - c], [d + e, +b]],
                 dtype=np.float64)
    E = np.array([[0, 1], [1, 2], [3, 4], [4, 5], [6, 7], [8, 9], [1, 4],
                  [3, 6], [5, 9]],
                 dtype=np.int32)

    return V, E


def main(args=None):
    parser = argparse.ArgumentParser(args)
    parser.add_argument("output", type=str, help="path to the output folder")
    parser.add_argument("--num-links", "-n", type=int, default=2)
    parser.add_argument("--iterations", "-m", type=int, default=20)
    args = parser.parse_args()

    output = args.output
    a = 0.4
    b = 0.5
    c = 0.40
    d = 0.4
    e = 0.5
    f = 0.25

    V, E = link(a, b, c, d, e)
    dx = np.array([d + e - f, 0.0], dtype=np.float64)
    rigid_bodies = []
    for i in range(0, args.num_links):
        dof_fixed = [False, False, False] if i != 0 else [True, True, True]

        rigid_bodies.append(
            dict(vertices=(V + dx * i).tolist(),
                 edges=E.tolist(),
                 velocity=[0.0, 0.0, 0.0],
                 is_dof_fixed=dof_fixed))
    data = {
        "max_iterations": args.iterations,
        "timestep": 0.1,
        "scene_type": "rigid_body_problem",
        "distance_barrier_constraint": {
            "initial_barrier_activation_distance": 0.1,
            "detection_method": "hash_grid",
            "use_distance_hashgrid": True,
            "custom_hashgrid_cellsize": -1
        },
        "barrier_solver": {
            "inner_solver": "newton_solver"
        },
        "rigid_body_problem": {
            "coefficient_restitution": 1.0,
            "gravity": [0.0, -0.5, 0.0],
            "rigid_bodies": rigid_bodies
        }
    }

    save_fixture(fixture, args.output)


if __name__ == "__main__":
    main()
