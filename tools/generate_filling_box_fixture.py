#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib

import numpy
import shapely.geometry
import shapely.ops

from fixture_utils import *


def generate_fixture(args):
    """Generate a fixture of a N boxes stacked on top of each other."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    hx = 5
    hy = 5
    half_thickness = 0.5
    box_polygons = [
        generate_rectangle(hx, half_thickness, numpy.array([0, 0]), 0),
        generate_rectangle(
            half_thickness, hy,
            numpy.array([hx - half_thickness, hy + half_thickness]), 0),
        generate_rectangle(
            half_thickness, hy,
            numpy.array([-hx + half_thickness, hy + half_thickness]), 0),
        generate_rectangle(
            hx, half_thickness,
            [hx - 2 * half_thickness, 2 * hy + half_thickness] +
            create_2D_rotation_matrix(numpy.pi / 4) @ [hx, -half_thickness],
            numpy.pi / 4),
        generate_rectangle(
            hx, half_thickness,
            [-hx + 2 * half_thickness, 2 * hy + half_thickness] +
            create_2D_rotation_matrix(-numpy.pi / 4) @ [-hx, -half_thickness],
            -numpy.pi / 4)
    ]

    box_vertices = list(
        shapely.ops.cascaded_union(box_polygons).exterior.coords)[:-1]
    box_polygons = [
        list(polygon.exterior.coords)[:-1] for polygon in box_polygons
    ]
    box_edges = generate_ngon_edges(len(box_vertices)).tolist()

    # Add the box
    rigid_bodies.append({
        "vertices": box_vertices,
        "polygons": box_polygons,
        "edges": box_edges,
        "oriented": False,
        "is_dof_fixed": [True, True, True],
    })

    radius = 1 / numpy.sqrt(2)
    block_vertices = [[0.5, -0.5], [0.5, 0.5], [-0.5, 0.5], [-0.5, -0.5]]

    centers = numpy.zeros((args.num_blocks, 2))

    # Add the pyramid
    for i in range(args.num_blocks):
        while True:
            center = numpy.random.random(2)
            centers - center
        rigid_bodies.append({
            "vertices": box_vertices,
            "polygons": [box_vertices],
            "edges": box_edges,
            "oriented": True,
            "position": center.tolist(),
            "theta": numpy.random.random() * numpy.pi / 4,
            "velocity": [0.0, 0.0, 0.0],
            "is_dof_fixed": [False, False, False],
        })

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser(description="generate a tower of blocks",
                                    default_initial_epsilon=1e-2,
                                    default_minimum_epsilon=1e-4,
                                    default_gravity=[0, -9.81, 0])
    parser.add_argument("--num-blocks",
                        type=int,
                        default=2,
                        help="number of blocks in the tower")
    args = parser.parse_args()

    if args.out_path is None:
        directory = pathlib.Path(__file__).resolve().parents[1] / "fixtures"
        args.out_path = (directory /
                         f"filling-box-num_blocks={args.num_blocks}.json")
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    with open(args.out_path, "w") as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
