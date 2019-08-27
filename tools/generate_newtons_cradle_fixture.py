#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import json
import numpy
import pathlib
import argparse

from default_fixture import generate_default_fixture


def generate_fixture(num_balls: int,
                     num_points: int,
                     walled: bool,
                     rotated: bool,
                     cor: float) -> dict:
    """Generate a fixture of a N boxes stacked on top of each other."""
    fixture = generate_default_fixture()
    fixture["distance_barrier_constraint"]["custom_initial_epsilon"] = 1e-2
    fixture["barrier_solver"]["min_barrier_epsilon"] = 1e-4
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    radius = 0.5
    x = numpy.cos(numpy.arange(num_points, dtype=float) /
                  num_points * 2 * numpy.pi) * radius
    y = numpy.sin(numpy.arange(num_points, dtype=float) /
                  num_points * 2 * numpy.pi) * radius
    vertices = numpy.hstack([x.reshape(-1, 1), y.reshape(-1, 1)]).tolist()
    edges = numpy.hstack([numpy.arange(num_points).reshape(-1, 1),
                          numpy.roll(numpy.arange(num_points).reshape(-1, 1), -1)]).tolist()

    if(walled):
        wall = {
            "vertices": [[0, -10], [0, 10]],
            "edges": [[0, 1]],
            "oriented": False,
            "theta": 0,
            "velocity": [0.0, 0.0, 0.0],
            "is_dof_fixed": [True, True, True]
        }
        wall["position"] = [-4.1, 0.0]
        rigid_bodies.append(wall.copy())
        wall["position"] = [num_balls * (2 * radius + 0.1), 0.0]
        rigid_bodies.append(wall)

    # Add the pyramid
    for i in range(num_balls):
        rigid_bodies.append({
            "vertices": vertices,
            "edges": edges,
            "oriented": True,
            "position": [-3 if i == 0 else ((i - 1) * (2 * radius + 0.1)), 0.0],
            "theta": ((i % 2) if rotated else 1) * numpy.rad2deg(numpy.pi / num_points),
            "velocity": [10.0 if i == 0 else 0.0, 0.0, 0.0],
            "is_dof_fixed": [False, False, False]
        })

    fixture["rigid_body_problem"]["gravity"] = [0, 0, 0]
    fixture["rigid_body_problem"]["coefficient_restitution"] = cor

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = argparse.ArgumentParser(
        description="generate a Newton's Cradle fixture")
    parser.add_argument("--num-balls", type=int, default=2,
                        help="number of balls in the cradle")
    parser.add_argument("--num-points", type=int, default=4,
                        help="number of points/edges used to discritize the balls")
    parser.add_argument("--not-walled", action="store_false", dest="walled",
                        help="do not put walls on both sides")
    parser.add_argument("--not-rotated", action="store_false", dest="rotated",
                        help="do not rotate alternating balls")
    parser.add_argument("--cor", type=float, default=-1,
                        help="coefficient of restitution")
    parser.add_argument("--out-path", metavar="path/to/output.json",
                        type=pathlib.Path, default=None,
                        help="path to save the fixture")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] /
                     "fixtures" / "newtons-cradle")
        directory.mkdir(parents=True, exist_ok=True)
        args.out_path = (directory /
                         "newtons-cradle-num_balls={:d}-num_points={:d}{}{}-cor={:g}.json".format(
                             args.num_balls, args.num_points,
                             "" if args.walled else "-not-walled",
                             "" if args.rotated else "-not-rotated", args.cor))
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print(args)

    fixture = generate_fixture(
        args.num_balls, args.num_points, args.walled, args.rotated, args.cor)

    with open(args.out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
