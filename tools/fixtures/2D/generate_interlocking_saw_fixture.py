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

    #  Add the saw
    num_teeth = args.num_teeth
    saw_length = 10
    tooth_length = saw_length / num_teeth

    saw_polygons = [
        (numpy.array([[0.0, 0.0], [tooth_length, 0.0], [tooth_length / 2, 1.0]
                      ]) + [tooth_length * i, 0]).tolist()
        for i in range(num_teeth)
    ]
    saw_vertices = numpy.array(saw_polygons).reshape(-1, 2).tolist()
    saw_edges = numpy.array([
        generate_ngon_edges(3) + 3 * i for i in range(num_teeth)
    ]).reshape(-1, 2).tolist()

    mass = 10  # Kg

    rigid_bodies.append({
        "vertices":
        saw_vertices,
        "polygons":
        saw_polygons,
        "edges":
        saw_edges,
        "oriented":
        True,
        "position": [0, 2],
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True],
        "masses":
        numpy.full(len(saw_vertices), mass / len(saw_vertices)).tolist(),
    })

    saw_polygons = [(numpy.array([[0.0, 0.0], [tooth_length / 2, -1.0],
                                  [tooth_length, 0.0]]) +
                     [tooth_length * (i + 0.5), 0]).tolist()
                    for i in range(num_teeth)]
    saw_vertices = numpy.array(saw_polygons).reshape(-1, 2).tolist()

    rigid_bodies.append({
        "vertices":
        saw_vertices,
        "polygons":
        saw_polygons,
        "edges":
        saw_edges,
        "oriented":
        True,
        "position": [0, 5],
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [False, False, False],
        "masses":
        numpy.full(len(saw_vertices), mass / len(saw_vertices)).tolist(),
    })

    # rigid_bodies.append(generate_walls_body(7, 7, [5, 4.5], 0.1))

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser("generate a block falling onto a saw",
                                    default_gravity=[0, -9.81, 0])
    parser.add_argument("--num-teeth",
                        type=int,
                        default=10,
                        help="number of teeth over a fixed length")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] / "fixtures" /
                     "saw")
        args.out_path = (
            directory /
            "interlocking-saws-num-teeth={:d}.json".format(args.num_teeth))
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
