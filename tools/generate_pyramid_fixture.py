#!/usr/local/bin/python3
"""Script to generate a fixture of a box falling on a saw."""

import sys
import json
import numpy
import pathlib

from default_fixture import generate_default_fixture


def generate_fixture():
    """Generate a fixture of a chain with N simple links."""
    fixture = generate_default_fixture()
    fixture["distance_barrier_constraint"]["custom_initial_epsilon"] = 1e-2
    fixture["barrier_solver"]["min_barrier_epsilon"] = 1e-4
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    # Add the ground
    rigid_bodies.append({
        "vertices": [[-10, 0], [10, 0]],
        "edges": [[0, 1]],
        "oriented": False,
        "is_dof_fixed": [True, True, True]
    })

    box_vertices = numpy.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    box_edges = [[0, 1], [1, 2], [2, 3], [3, 0]]

    # Add the pyramid
    nrows = 5
    for i in range(nrows):
        y_offset = (1 + 1e-1) * i + 1e-1
        for j in range(nrows - i):
            x_offset = (1 + 1e-1) * j - (1 + 1e-1) * (nrows - i) / 2
            rigid_bodies.append({
                "vertices": (box_vertices + [x_offset, y_offset]).tolist(),
                "edges": box_edges,
                "velocity": [0.0, 0.0, 0.0],
                "is_dof_fixed": [False, False, False]
            })

    fixture["rigid_body_problem"]["gravity"] = [0, -9.81, 0]
    fixture["rigid_body_problem"]["coefficient_restitution"] = -1

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    fixture = generate_fixture()

    if(len(sys.argv) > 1):
        out_path = pathlib.Path(sys.argv[2])
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        directory = (pathlib.Path(__file__).resolve().parents[1] /
                     "fixtures" / "stacking")
        directory.mkdir(parents=True, exist_ok=True)
        out_path = directory / "pyramid.json"
    with open(out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
