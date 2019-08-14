#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import sys
import json
import numpy
import pathlib

from default_fixture import generate_default_fixture


def generate_fixture(height=2):
    """Generate a fixture of a N boxes stacked on top of each other."""
    fixture = generate_default_fixture()
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    box_vertices = numpy.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    box_edges = [[0, 1], [1, 2], [2, 3], [3, 0]]

    # Add the ground
    rigid_bodies.append({
        "vertices": [[-10, -1], [10, -1], [10, 0], [-10, 0]],
        "edges": box_edges,
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    })

    # Add the pyramid
    for i in range(height):
        rigid_bodies.append({
            "vertices": (box_vertices + [0, (1 + 1e-1) * i + 1e-1]).tolist(),
            "edges": box_edges,
            "velocity": [0.0, 0.0, 0.0],
            "is_dof_fixed": [False, False, False]
        })

    fixture["rigid_body_problem"]["gravity"] = [0, -9.81, 0]

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    if len(sys.argv) > 1:
        height = int(sys.argv[1])
    else:
        height = 2

    fixture = generate_fixture(height)

    if len(sys.argv) > 2:
        out_path = pathlib.Path(sys.argv[2])
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        directory = (pathlib.Path(__file__).resolve().parents[1] /
                     "fixtures" / "rigid_bodies")
        directory.mkdir(parents=True, exist_ok=True)
        out_path = directory / f"stacking-{height}-tower.json"

    with open(out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
