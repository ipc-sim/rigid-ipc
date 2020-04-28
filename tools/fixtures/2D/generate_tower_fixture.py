#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib

import sys

import numpy

from fixture_utils import *


def generate_fixture(args):
    """Generate a fixture of a N boxes stacked on top of each other."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    # Generate tower
    radius = 1 / numpy.sqrt(2)
    box_vertices = [[0.5, -0.5], [0.5, 0.5], [-0.5, 0.5], [-0.5, -0.5]]
    box_edges = [[0, 1], [1, 2], [2, 3], [3, 0]]
    for i in range(args.num_blocks):
        x = args.x_offset if i % 2 else 0
        y = (2 if i == args.num_blocks - 1 and args.falling else
             1) * 2 * radius * i + 1.25 * radius
        theta = 0 if not args.rotated or i % 2 else 45
        rigid_bodies.append({
            "vertices": box_vertices,
            "polygons": [box_vertices],
            "edges": box_edges,
            "oriented": True,
            "position": [x, y],
            "theta": theta,
            "velocity": [0.0, 0.0, 0.0],
            "is_dof_fixed": [False, False, False],
            "masses": numpy.full(4, 0.25).tolist(),
            "density": 1
        })

    # Add the walls around tower
    y += 2
    rigid_bodies.append(
        generate_walls_body(5, y / 2, numpy.array([0, y / 2]), 0.1))

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    print(sys.argv)
    parser = create_argument_parser(description="generate a tower of blocks",
                                    default_initial_epsilon=1e-2,
                                    default_gravity=[0, -9.81, 0])
    parser.add_argument("--num-blocks",
                        type=int,
                        default=2,
                        help="number of blocks in the tower")
    parser.add_argument("--x-offset",
                        type=float,
                        default=0,
                        help="offset alternating blocks in x")
    parser.add_argument("--rotated",
                        action="store_true",
                        help="rotate alternating blocks")
    parser.add_argument("--falling",
                        action="store_true",
                        help="last block falling from high")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] / "fixtures" /
                     "stacking")
        args.out_path = (
            directory /
            "tower-num_blocks={:d}-cor={:g}-x_offset={:g}{}{}.json".format(
                args.num_blocks, args.restitution_coeff, args.x_offset,
                "-rotated" if args.rotated else "",
                "-falling" if args.falling else ""))
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
