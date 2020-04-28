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

    # Add the ground
    rigid_bodies.append(generate_walls_body(10, 4, [0, 4], 0.1))

    # Add the pyramid
    box = generate_box_body(0.5, 0.5, [0, 0], 0, 10)
    nrows = 5
    for i in range(nrows):
        y = (1 + 1e-1) * i + 0.6
        for j in range(nrows - i):
            x = (1 + 1e-1) * j - (1 + 1e-1) * (nrows - i) / 2
            box["position"] = [x, y]
            rigid_bodies.append(box.copy())

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

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
