#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import pathlib

import numpy

from fixture_utils import *


def generate_fixture(args):
    """Generate a fixture of a N boxes stacked on top of each other."""
    fixture = generate_custom_fixture(args)
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    # Generate lines
    line_vertices = generate_rectangle_vertices(5, args.thickness / 2, [0, 0],
                                                0)
    line_edges = generate_ngon_edges(4)
    line = {
        "vertices": line_vertices.tolist(),
        "polygons": [line_vertices.tolist()],
        "edges": line_edges.tolist(),
        "oriented": True,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [False, False, False]
    }
    delta_y = 0.5 + args.thickness
    y = -args.thickness / 2
    for i in range(args.num_lines):
        y += delta_y
        line["position"] = [0, y]
        rigid_bodies.append(line.copy())

    # Add box falling
    y += 10
    line["vertices"] = (line_vertices * [0.1, 1 / args.thickness]).tolist()
    line["polygons"] = [line["vertices"]]
    line["position"] = [0, y]
    line["theta"] = 10  # °
    line["velocity"] = [0, args.impact_velocity, 0]
    rigid_bodies.append(line.copy())

    # Add the walls around tower
    y += 2
    rigid_bodies.append(
        generate_walls(numpy.array([0, y / 2]), 5.5, y / 2, 0.1))

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    parser = create_argument_parser(description="generate a stack of lines",
                                    default_gravity=[0, 0, 0])
    parser.add_argument("--num-lines",
                        type=int,
                        default=10,
                        help="number of lines in the stack")
    parser.add_argument("--thickness",
                        type=float,
                        default=1e-2,
                        help="thickness of each line")
    parser.add_argument("--impact-velocity",
                        type=float,
                        default=-100,
                        help="thickness of each line")
    args = parser.parse_args()

    if args.out_path is None:
        directory = (pathlib.Path(__file__).resolve().parents[1] / "fixtures" /
                     "stacking")
        args.out_path = (
            directory /
            "line-stack-num_line={:d}-thickness={:g}-v0={:g}.json".format(
                args.num_lines, args.thickness, args.impact_velocity))
    args.out_path.parent.mkdir(parents=True, exist_ok=True)

    print_args(args)

    fixture = generate_fixture(args)

    save_fixture(fixture, args.out_path)


if __name__ == "__main__":
    main()
