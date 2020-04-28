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
        "theta": 22.5,
        "is_dof_fixed": [True, True, True]
    })

    num_points = 25  # number of points on the interior/exterior

    # Set the radii of the wheel
    wheel_inner_radius = axle_radius + 0.5
    wheel_outer_radius = wheel_inner_radius + 1

    # Vertices of the wheel
    wheel_inner_vertices = generate_regular_ngon_vertices(
        num_points, wheel_inner_radius)[::-1]  # Inner verties should be CW
    wheel_outer_vertices = generate_regular_ngon_vertices(
        num_points, wheel_outer_radius)  # Outer verties should be CW
    wheel_vertices = numpy.vstack([wheel_inner_vertices, wheel_outer_vertices])
    # Edges of the wheel
    wheel_edges = generate_ngon_edges(num_points)
    wheel_edges = numpy.vstack([wheel_edges, wheel_edges + num_points])
    # Decompose the wheel into convex quadralaterals
    wheel_polygons = numpy.array([[
        wheel_inner_vertices[i], wheel_inner_vertices[(i + 1) % num_points],
        wheel_outer_vertices[(num_points - (i + 2) % num_points) % num_points],
        wheel_outer_vertices[(num_points - (i + 1) % num_points) % num_points]
    ] for i in range(num_points)])
    for polygon in wheel_polygons:
        assert is_polygon_ccw(polygon)

    mass = 100  # Kg
    area = (compute_regular_ngon_area(wheel_outer_vertices) -
            compute_regular_ngon_area(wheel_inner_vertices))  # m²
    density = mass / area  # Kg / m²

    rigid_bodies.append({
        "vertices":
        wheel_vertices.tolist(),
        "polygons":
        wheel_polygons.tolist(),
        "edges":
        wheel_edges.tolist(),
        "oriented":
        True,
        "velocity": [0.0, 0.0, 10 * numpy.pi],
        "density":
        density,
        "masses":
        numpy.full(wheel_vertices.shape[0],
                   mass / wheel_vertices.shape[0]).tolist(),
        "is_dof_fixed": [False, False, False]
    })

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser(
        "generate a wheel spinning loose on an axle",
        default_gravity=[0, -9.81, 0])
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
