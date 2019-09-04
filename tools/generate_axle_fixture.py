#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib
import sys

import numpy
import shapely.geometry
import shapely.ops

from fixture_utils import *


def generate_fixture(args):
    """Generate a saw and block."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    axle_radius = 0.5
    axle_vertices = generate_regular_ngon_vertices(8, axle_radius)
    axle_edges = generate_ngon_edges(8)

    rigid_bodies.append({
        "vertices": axle_vertices.tolist(),
        "polygons": [axle_vertices.tolist()],
        "edges": axle_edges.tolist(),
        "oriented": True,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    })

    num_points = 8
    wheel_inner_radius = axle_radius + 0.5
    wheel_outer_radius = annulus_inner_radius + 1
    scaling = wheel_outer_radius / wheel_inner_radius
    wheel_inner_vertices = generate_regular_ngon_vertices(
        num_points, wheel_inner_radius)
    wheel_outer_radius = scaling * wheel_inner_radius

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser("generate a block falling onto a saw",
                                    default_minimum_epsilon=1e-4,
                                    default_gravity=[0, -9.8, 0])
    args = parser.parse_args()

    if args.out_path is None:
        directory = pathlib.Path(__file__).resolve().parents[1] / "fixtures"
        args.out_path = directory / "axle.json"
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
