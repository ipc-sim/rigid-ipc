#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib

import numpy

from fixture_utils import *


def generate_fixture(args):
    """Generate a fixture of a chain with N simple links."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    box_edges = generate_ngon_edges(4).tolist()

    # Add the ground
    ground_box_vertices = [[10, 0], [-10, 0], [-10, -1], [10, -1]]
    rigid_bodies.append({
        "vertices": ground_box_vertices,
        "polygons": [ground_box_vertices],
        "edges": box_edges,
        "oriented": False,
        "is_dof_fixed": [True, True, True]
    })

    box_vertices = [[0, 0], [1, 0], [1, 1], [0, 1]]

    # Add the pyramid
    nrows = 5
    for i in range(nrows):
        y = (1 + 1e-1) * i + 1e-1
        for j in range(nrows - i):
            x = (1 + 1e-1) * j - (1 + 1e-1) * (nrows - i) / 2
            rigid_bodies.append({
                "vertices": box_vertices,
                "polygons": [box_vertices],
                "edges": box_edges,
                "position": [x, y],
                "velocity": [0.0, 0.0, 0.0],
                "is_dof_fixed": [False, False, False]
            })

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser("generate a pyramid fixture",
                                    default_gravity=[0, -9.81, 0])
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] / "fixtures" /
                     "stacking")
        args.out_path = (
            directory / "pyramid-cor={:g}.json".format(args.restitution_coeff))
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    with open(args.out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
