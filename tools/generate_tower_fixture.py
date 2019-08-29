#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import numpy
import pathlib

from default_fixture import generate_default_fixture


def generate_fixture(num_blocks, cor, x_offset, rotated, falling):
    """Generate a fixture of a N boxes stacked on top of each other."""
    fixture = generate_default_fixture()
    fixture["distance_barrier_constraint"]["custom_initial_epsilon"] = 1e-2
    fixture["barrier_solver"]["min_barrier_epsilon"] = 1e-4
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    # Add the ground
    rigid_bodies.append({
        "vertices": [[-10, 0], [10, 0]],
        "edges": [[0, 1]],
        "oriented": False,
        "is_dof_fixed": [True, True, True]
    })

    radius = 1 / numpy.sqrt(2)
    box_vertices = [[0.5, -0.5], [0.5, 0.5], [-0.5, 0.5], [-0.5, -0.5]]
    box_edges = [[0, 1], [1, 2], [2, 3], [3, 0]]

    # Add the pyramid
    for i in range(num_blocks):
        x = x_offset if i % 2 else 0
        y = ((2 if i == num_blocks - 1 and falling else 1)
             * 2 * radius * i + 1.25 * radius)
        theta = 0 if not rotated or i % 2 else 45
        rigid_bodies.append({
            "vertices": box_vertices,
            "edges": box_edges,
            "oriented": True,
            "position": [x, y],
            "theta": theta,
            "velocity": [0.0, 0.0, 0.0],
            "is_dof_fixed": [False, False, False]
        })

    fixture["rigid_body_problem"]["gravity"] = [0, -9.81, 0]
    fixture["rigid_body_problem"]["coefficient_restitution"] = cor

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = argparse.ArgumentParser(
        description="generate a tower of blocks")
    parser.add_argument("--num-blocks", type=int, default=2,
                        help="number of blocks in the tower")
    parser.add_argument("--cor", type=float, default=-1,
                        help="coefficient of restitution")
    parser.add_argument("--x-offset", type=float, default=0,
                        help="offset alternating blocks in x")
    parser.add_argument("--rotated", action="store_true",
                        help="rotate alternating blocks")
    parser.add_argument("--falling", action="store_true",
                        help="last block falling from high")
    parser.add_argument("--out-path", metavar="path/to/output.json",
                        type=pathlib.Path, default=None,
                        help="path to save the fixture")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] /
                     "fixtures" / "stacking")
        args.out_path = (
            directory / "tower-num_blocks={:d}-cor={:g}-x_offset={:g}{}{}.json".format(
                args.num_blocks, args.cor, args.x_offset,
                "-rotated" if args.rotated else "",
                "-falling" if args.falling else ""))
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print(args)

    fixture = generate_fixture(
        args.num_blocks, args.cor, args.x_offset, args.rotated, args.falling)

    with open(args.out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
