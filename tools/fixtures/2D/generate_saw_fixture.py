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

    # Add the box
    box = generate_box_body(0.5, 0.5, [-0.5, 2.5], 0, 10)
    box["velocity"] = [1, -1, 0]
    rigid_bodies.append(box)

    #  Add the saw
    num_teeth = 100
    saw_length = 10
    tooth_length = saw_length / num_teeth
    saw_polygons = [
        shapely.geometry.Polygon(
            numpy.array([[0.0, 0.0], [tooth_length, 0.0], [0.0, 1.0]]) +
            [0.8 * tooth_length * i, 0]) for i in range(num_teeth)
    ]
    saw = shapely.ops.cascaded_union(saw_polygons)
    saw = shapely.geometry.polygon.orient(saw, 1)
    saw_polygons = [list(polygon.exterior.coords) for polygon in saw_polygons]
    saw_vertices = numpy.array(list(saw.exterior.coords)[:-1])
    saw_edges = generate_ngon_edges(saw_vertices.shape[0])

    rigid_bodies.append({
        "vertices": saw_vertices.tolist(),
        "polygons": saw_polygons,
        "edges": saw_edges.tolist(),
        "oriented": True,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    })

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser("generate a block falling onto a saw")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] / "fixtures" /
                     "saw")
        args.out_path = (directory /
                         "saw-cor={:g}.json".format(args.restitution_coeff))
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
