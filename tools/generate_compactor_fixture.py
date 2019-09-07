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
    rigid_bodies.append(generate_walls_body(5, 2.5, numpy.zeros(2), 0.1))

    # Compactor Box
    hx = hy = 2.49
    compactor_mass = 1000  # kg
    compactor = generate_box_body(hx, hy, [-hx, 0], 0, compactor_mass)
    compactor["velocity"][0] = 100
    rigid_bodies.append(compactor)

    # Trash Boxes
    radius = 0.25
    hx = hy = numpy.sqrt(radius**2 / 2)
    radius += 5e-2  # inflate the radius slightly
    trash = generate_box_body(hx, hy, [0, 0], 0, 1)

    num_trash = args.num_blocks
    centers = numpy.zeros((num_trash, 2))

    width = 5 - 2 * radius
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
                      [radius, -height / 2])
            invalid_center = (numpy.linalg.norm(centers - center, axis=1) <
                              (2 * radius)).any()
            num_tries += 1

        centers[i] = center
        trash["position"] = center.tolist()
        trash["theta"] = numpy.random.random() * 45
        rigid_bodies.append(trash.copy())

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    print(sys.argv)
    parser = create_argument_parser(description="generate a tower of blocks",
                                    default_restitution_coefficient=0)
    parser.add_argument("--num-blocks",
                        type=int,
                        default=10,
                        help="maximum number of blocks in the compactor")
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
