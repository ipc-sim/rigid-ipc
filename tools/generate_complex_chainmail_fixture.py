#!/usr/local/bin/python3
"""
Script to generate a fixture of a chain with N complex links.

Usage: python generate_complex_chainmail_fixture.py N
"""

import sys
import json
import numpy
import pathlib

from default_fixture import generate_default_fixture


def generate_fixture(n_links: int) -> dict:
    """Generate a fixture of a chain with N complex links."""
    link_geometry_path = (pathlib.Path(__file__).resolve().parent /
                          "complex-link.json")
    with open(link_geometry_path) as link_file:
        link_mesh = json.load(link_file)
    vertices = numpy.array(link_mesh["vertices"], dtype=float)
    edges = numpy.array(link_mesh["edges"], dtype=int)
    fixture = generate_default_fixture()
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]
    for i in range(n_links):
        rigid_bodies.append({
            "vertices": (vertices + [0, -1.9 * i]).tolist(),
            "edges": edges.tolist(),
            "velocity": [0.0, -1.0 if i else 0.0, 0.0],
            "is_dof_fixed": numpy.full((3,), i == 0, dtype=bool).tolist()
        })
    return fixture


def main() -> None:
    """Parse command-line arguments to generate the desired fixture."""
    assert(len(sys.argv) >= 2)

    n_links = int(sys.argv[1])
    fixture = generate_fixture(n_links)

    if(len(sys.argv) > 2):
        out_path = pathlib.Path(sys.argv[2])
        out_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        directory = (pathlib.Path(__file__).resolve().parents[1] /
                     "fixtures" / "chain")
        directory.mkdir(parents=True, exist_ok=True)
        out_path = directory / f"complex_{n_links:d}_link_chain.json"
    with open(out_path, 'w') as outfile:
        json.dump(fixture, outfile)


if __name__ == "__main__":
    main()
