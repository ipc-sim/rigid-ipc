#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib

import numpy
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

from fixture_utils import *


def generate_fixture(args: argparse.Namespace) -> dict:
    """Generate a fixture of a N boxes stacked on top of each other."""
    numpy.random.seed(seed=0)  # Deterministic random results

    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    # Dimensions of the table:
    # https://www.dimensions.guide/element/7-foot-billiards-pool-table
    width = 1.98  # m
    height = .99  # m
    thickness = .0019  # m
    rigid_bodies.append(
        generate_walls_body(width / 2, height / 2, numpy.zeros(2), thickness))

    # Ball settings
    radius = 0.057 / 2  # m
    ball_vertices = generate_regular_ngon_vertices(args.num_points, radius)

    mass = 0.170  # Kg
    area = compute_regular_ngon_area(ball_vertices)  # m²
    density = mass / area  # Kg / m²

    ball = {
        "vertices": ball_vertices.tolist(),
        "polygons": [ball_vertices.tolist()],
        "edges": generate_ngon_edges(args.num_points).tolist(),
        "oriented": True,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [False, False, False],
        "masses": numpy.full(args.num_points, mass / args.num_points).tolist(),
        "density": density
    }

    # Add the triangle of balls
    x = width / 4
    for row in range(5):
        y = -1.15 * radius * row
        for i in range(row + 1):
            ball["position"] = [x, y]
            ball["theta"] = numpy.random.random() * (360 / args.num_points)
            rigid_bodies.append(ball.copy())
            y += 2.3 * radius
        x += 2 * radius

    # Add the cue ball
    ball["position"] = [-width / 4, 0]
    ball["theta"] = 0
    # Typical Pool Ball Speeds: Powerful Break
    # https://billiards.colostate.edu/faq/speed/typical/
    break_speed = 10.729  # m / s
    ball["velocity"] = [break_speed, 0.0, 0.0]
    rigid_bodies.append(ball.copy())

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser("generate a billiards fixture",
                                    default_timestep=1e-3,
                                    default_initial_epsilon=1e-2,
                                    default_restitution_coefficient=1)
    parser.add_argument(
        "--num-points",
        type=int,
        default=8,
        help="number of points/edges used to discritize the balls")
    args = parser.parse_args()

    if args.out_path is None:
        directory = pathlib.Path(__file__).resolve().parents[1] / "fixtures"
        args.out_path = directory / "billiards.json"
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
