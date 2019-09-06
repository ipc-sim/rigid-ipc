#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import argparse
import json
import pathlib

import sys

import numpy

from fixture_utils import *


def generate_fixture(args):
    """Generate a fixture of a N boxes stacked on top of each other."""
    numpy.random.seed(0)
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    # Add the walls around the scene
    rigid_bodies.append(generate_walls(numpy.array([0, 0]), 5, 2.5, 0.1))

    # Compactor Box
    hx = hy = 2.49
    compactor_vertices = generate_rectangle_vertices(hx, hy, [0, 0], 0)
    compactor_edges = generate_ngon_edges(4)
    compactor_mass = 1000  # kg
    compactor_area = 4 * hx * hy  # m²
    compactor_density = compactor_mass / compactor_area  # kg / m²
    rigid_bodies.append({
        "vertices": compactor_vertices.tolist(),
        "polygons": [compactor_vertices.tolist()],
        "edges": compactor_edges.tolist(),
        "oriented": True,
        "position": [-hx, 0],
        "velocity": [100.0, 0.0, 0.0],
        "masses": numpy.full(4, compactor_mass / 4).tolist(),
        "density": compactor_density,
        "is_dof_fixed": [False, False, False]
    })

    # Trash Boxes
    radius = 0.25
    trash_vertices = generate_regular_ngon_vertices(4, radius).tolist()
    trash_edges = generate_ngon_edges(4).tolist()

    num_trash = 10
    centers = numpy.zeros((num_trash, 2))

    width = 5 - 2 * radius - 1e-1
    height = width
    for i in range(num_trash):
        invalid_center = True
        num_tries = 0
        while invalid_center:
            if num_tries > 1000:
                print("Less than 0.1% chance of finding new center.")
                print(f"Quiting with a maximum of {i:d} trash blocks.")
                return fixture
            center = (numpy.random.random(2) * [width, height] +
                      [radius + 5e-2, -height / 2])
            invalid_center = (numpy.linalg.norm(centers - center, axis=1) <
                              (2 * radius + 1e-4)).any()
            num_tries += 1

        centers[i] = center
        rigid_bodies.append({
            "vertices": trash_vertices,
            "polygons": [trash_vertices],
            "edges": trash_edges,
            "oriented": True,
            "position": center.tolist(),
            "theta": numpy.random.random() * 45,
            "velocity": [0.0, 0.0, 0.0],
            "is_dof_fixed": [False, False, False],
        })

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    print(sys.argv)
    parser = create_argument_parser(description="generate a tower of blocks",
                                    default_restitution_coefficient=0)
    args = parser.parse_args()

    if args.out_path is None:
        directory = pathlib.Path(__file__).resolve().parents[1] / "fixtures"
        args.out_path = directory / "compactor.json"
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
