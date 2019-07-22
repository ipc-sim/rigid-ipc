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
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]

    # Add the box
    rigid_bodies.append({
        "vertices": [[2, 3], [0, 3], [0, 2], [2, 2]],
        "edges": numpy.vstack(
            [numpy.arange(4), numpy.roll(numpy.arange(4), -1)]).T.tolist(),
        "velocity": [1.0, -1.0, 0.0]
    })

    #  Add the saw
    vertices = []
    min_x, max_x, delta_x = 1.0, 4.0, 0.1
    x = min_x
    while x <= max_x:
        vertices.append([x, 1.0])
        x = round(x + delta_x, 4)
        vertices.append([x, 0.0])
    rigid_bodies.append({
        "vertices": vertices,
        "edges": numpy.vstack([
            numpy.arange(len(vertices)),
            numpy.roll(numpy.arange(len(vertices)), -1)]).T[:-1].tolist(),
        "velocity": [0.0, 0.0, 0.0],
        "is_dof_fixed": [True, True, True]
    })

    return fixture


def main():
    """Parse command-line arguments to generate the desired fixture."""
    fixture = generate_fixture()

    if(len(sys.argv) > 1):
        out_path = pathlib.Path(sys.argv[2])
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        directory = (pathlib.Path(__file__).resolve().parents[1] /
                     "fixtures" / "rigid_bodies")
        directory.mkdir(parents=True, exist_ok=True)
        out_path = directory / "square_falling_onto_saw.json"
    with open(out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
