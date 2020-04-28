#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib

import numpy

from fixture_utils import *


def generate_fixture(args: argparse.Namespace) -> dict:
    """Generate a fixture of Newton's cradle."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    radius = 0.5

    hx = (args.num_balls * (2 * radius + 0.1) + 4.1) / 2
    cx = hx - 4.1
    rigid_bodies.append(generate_walls_body(hx, hx, numpy.array([cx, 0]), 0.1))

    ball_vertices = generate_regular_ngon_vertices(args.num_points, radius)
    ball_edges = generate_ngon_edges(args.num_points)

    ball_mass = 1
    ball_area = compute_regular_ngon_area(ball_vertices)
    ball_density = ball_mass / ball_area

    ball = {
        "vertices": ball_vertices.tolist(),
        "polygons": [ball_vertices.tolist()],
        "edges": ball_edges.tolist(),
        "oriented": True,
        "is_dof_fixed": [False, False, False],
        "masses": numpy.full(args.num_points,
                             ball_mass / args.num_points).tolist(),
        "density": ball_density
    }
    for i in range(args.num_balls):
        ball["position"] = [
            -3 if i == 0 else ((i - 1) * (2 * radius + 0.1)), 0.0
        ]
        ball["theta"] = ((i % 2) if args.rotated else 1) * numpy.rad2deg(
            numpy.pi / args.num_points)
        ball["velocity"] = [10.0 if i == 0 else 0.0, 0, 0]
        rigid_bodies.append(ball.copy())

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser(
        description="generate a Newton's Cradle fixture",
        default_initial_epsilon=1e-2,
        default_restitution_coefficient=1)
    parser.add_argument("--num-balls",
                        type=int,
                        default=5,
                        help="number of balls in the cradle")
    parser.add_argument(
        "--num-points",
        type=int,
        default=8,
        help="number of points/edges used to discritize the balls")
    parser.add_argument("--not-rotated",
                        action="store_false",
                        dest="rotated",
                        help="do not rotate alternating balls")
   
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] / "fixtures" /
                     "newtons-cradle")
        filename = ("newtons-cradle"
                    f"-num_balls={args.num_balls:d}"
                    f"-num_points={args.num_points:d}"
                    f"{'' if args.rotated else '-not-rotated':s}"
                    f"-cor={args.restitution_coeff:g}"
                    ".json")
        args.out_path = directory / filename
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
