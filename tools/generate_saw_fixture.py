#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib
import sys

import numpy
import shapely.geometry
import shapely.ops

from default_fixture import generate_default_fixture


def generate_fixture(cor):
    """Generate a saw and block."""
    fixture = generate_default_fixture()
    fixture["distance_barrier_constraint"]["custom_initial_epsilon"] = 1e-1
    fixture["barrier_solver"]["min_barrier_epsilon"] = 1e-4
    fixture["rigid_body_problem"]["coefficient_restitution"] = cor
    fixture["rigid_body_problem"]["gravity"] = [0, 0, 0]
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    # Add the box
    rigid_bodies.append({
        "vertices": [[-1, 2], [0, 2], [0, 3], [-1, 3]],
        "polygons": [[[-1, 2], [0, 2], [0, 3], [-1, 3]]],
        "edges": numpy.vstack(
            [numpy.arange(4), numpy.roll(numpy.arange(4), -1)]).T.tolist(),
        "velocity": [1.0, -1.0, 0.0]
    })

    #  Add the saw
    num_teeth = 100
    saw_length = 10
    tooth_length = saw_length / num_teeth
    saw_polygons = [shapely.geometry.Polygon(numpy.array(
        [[0.0, 0.0], [tooth_length, 0.0], [0.0, 1.0]]) +
        [0.8 * tooth_length * i, 0])
        for i in range(num_teeth)]
    saw = shapely.ops.cascaded_union(saw_polygons)
    saw_polygons = [list(polygon.exterior.coords)
                    for polygon in saw_polygons]
    vertices = numpy.array(list(saw.exterior.coords)[:-1])
    num_points = vertices.shape[0]
    edges = numpy.hstack([
        numpy.arange(num_points).reshape(-1, 1),
        numpy.roll(numpy.arange(num_points).reshape(-1, 1), -1)])

    rigid_bodies.append({
        "vertices": vertices.tolist(),
        "edges": edges.tolist(),
        "polygons": saw_polygons,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    })

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = argparse.ArgumentParser(
        description="generate a block falling onto a saw")
    parser.add_argument("--cor", type=float, default=-1,
                        help="coefficient of restitution")
    parser.add_argument("--out-path", metavar="path/to/output.json",
                        type=pathlib.Path, default=None,
                        help="path to save the fixture")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] /
                     "fixtures" / "saw")
        args.out_path = directory / "saw-cor={:g}.json".format(args.cor)
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print(args)

    fixture = generate_fixture(args.cor)

    with open(args.out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
