#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib
import sys

import numpy

from fixture_utils import *


def generate_fixture(args):
    """Generate a saw and block."""
    fixture = generate_custom_fixture(args)
    mass = 1
    area = 4
    density = mass / area
    fixture["rigid_body_problem"]["rigid_bodies"] = [{
        "vertices": [[-1, -1], [1, -1], [1, 1], [-1, 1]],
        "polygons": [[[-1, -1], [1, -1], [1, 1], [-1, 1]]],
        "edges": [[0, 1], [1, 2], [2, 3], [3, 0]],
        "oriented":
        True,
        "position": [-10, 10],
        "theta":
        45,
        "velocity": [1.0, 0.0, 0.0],
        "is_dof_fixed": [False, False, False],
        "masses":
        numpy.full(4, mass / 4).tolist(),
        "density":
        density
    }, {
        "vertices": [[-10, 0], [10, 0], [10, 1], [-10, 1]],
        "polygons": [[[-10, 0], [10, 0], [10, 1], [-10, 1]]],
        "edges": [[0, 1], [1, 2], [2, 3], [3, 0]],
        "oriented":
        True,
        "position": [0, 0],
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    }]

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser(
        "generate a wheel spinning loose on an axle",
        default_gravity=[0, -3.0, 0],
        default_num_steps=2300,
        default_timestep=1e-2,
        default_restitution_coefficient=1)
    args = parser.parse_args()

    if args.out_path is None:
        directory = pathlib.Path(__file__).resolve().parents[1] / "fixtures"
        args.out_path = directory / "bouncing-diamond.json"
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
