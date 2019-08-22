#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import json
import numpy
import pathlib
import argparse

from default_fixture import generate_default_fixture


def generate_fixture(num_points: int):
    """Generate a fixture of a N boxes stacked on top of each other."""
    fixture = generate_default_fixture()
    fixture["timestep_size"] = 1e-3
    fixture["distance_barrier_constraint"]["custom_initial_epsilon"] = 1e-2
    fixture["barrier_solver"]["min_barrier_epsilon"] = 1e-4
    fixture["rigid_body_problem"]["gravity"] = [0, 0, 0]
    fixture["rigid_body_problem"]["coefficient_restitution"] = 1
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    width = 198  # cm
    height = 99  # cm
    # Walls of the pool table
    rigid_bodies.append({
        "vertices": [[-width / 2, -height / 2], [width / 2, -height / 2],
                     [width / 2, height / 2], [-width / 2, height / 2]],
        "edges": [[0, 1], [1, 2], [2, 3], [3, 0]],
        "oriented": False,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    })

    # Ball settings
    radius = 5.7 / 2  # cm
    x = numpy.cos(numpy.arange(num_points, dtype=float) /
                  num_points * 2 * numpy.pi) * radius
    y = numpy.sin(numpy.arange(num_points, dtype=float) /
                  num_points * 2 * numpy.pi) * radius
    ball = {
        "vertices": numpy.hstack([x.reshape(-1, 1), y.reshape(-1, 1)]).tolist(),
        "edges": numpy.hstack([numpy.arange(num_points).reshape(-1, 1),
                               numpy.roll(numpy.arange(num_points).reshape(-1, 1), -1)]).tolist(),
        "oriented": True,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [False, False, False]
    }

    # Add the triangle of balls
    x = width / 4
    for row in range(5):
        y = -1.15 * radius * row
        for i in range(row + 1):
            ball["position"] = [x, y]
            ball["theta"] = numpy.random.random() * (360 / num_points)
            rigid_bodies.append(ball.copy())
            y += 2.3 * radius
        x += 2 * radius

    # Add the cue ball
    ball["position"] = [-width / 4, 0]
    ball["theta"] = 0
    break_speed = 1072.9  # cm / s
    ball["velocity"] = [break_speed, 0.0, 0.0]
    rigid_bodies.append(ball.copy())

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = argparse.ArgumentParser(
        description="generate a Newton's Cradle fixture")
    parser.add_argument("--num-points", type=int, default=25,
                        help="number of points/edges used to discritize the balls")
    parser.add_argument("--out-path", metavar="path/to/output.json",
                        type=pathlib.Path, default=None,
                        help="path to save the fixture")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] /
                     "fixtures")
        directory.mkdir(parents=True, exist_ok=True)
        args.out_path = directory / "billiards.json"
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print(args)

    fixture = generate_fixture(args.num_points)

    with open(args.out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
