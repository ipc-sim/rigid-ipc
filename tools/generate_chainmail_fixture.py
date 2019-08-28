#!/usr/local/bin/python3
"""
Script to generate a fixture of a chain with N complex links.

Usage: python generate_complex_chainmail_fixture.py N
"""

import sys
import json
import numpy
import pathlib
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

from default_fixture import generate_default_fixture


def create_link_polygons():
    return [
        Polygon([(1.5, 6), (4.5, 6), (4.5, 7), (1.5, 7)]),
        Polygon([(2.5, 3.5), (3.5, 3.5), (3.5, 6.5), (2.5, 6.5)]),
        Polygon([(0, 3), (6, 3), (6, 4), (0, 4)]),
        Polygon([(0, 0.5), (1, 0.5), (1, 3.5), (0, 3.5)]),
        Polygon([(0, 0), (2, 0), (2, 1), (0, 1)]),
        Polygon([(5, 0.5), (6, 0.5), (6, 3.5), (5, 3.5)]),
        Polygon([(4, 0), (6, 0), (6, 1), (4, 1)]),
    ]


def generate_fixture(n_links: int) -> dict:
    """Generate a fixture of a chain with N complex links."""
    link_polygons = create_link_polygons()
    link = cascaded_union(link_polygons)
    link_polygons = [list(polygon.exterior.coords)
                     for polygon in link_polygons]
    vertices = numpy.array(list(link.exterior.coords)[:-1])

    angle = 90 + 45
    theta = numpy.radians(angle)
    c, s = numpy.cos(theta), numpy.sin(theta)
    R = numpy.array(((c, -s), (s, c)))

    num_points = vertices.shape[0]
    edges = numpy.hstack([
        numpy.arange(num_points).reshape(-1, 1),
        numpy.roll(numpy.arange(num_points).reshape(-1, 1), -1)]).tolist()

    fixture = generate_default_fixture()
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]
    for i in range(n_links):
        rigid_bodies.append({
            "vertices": vertices.tolist(),
            "edges": edges,
            "polygons": link_polygons,
            "position": (R @ numpy.array([0, 6 * n_links - 4.5 * i])).tolist(),
            "theta": angle,
            "velocity": [0.0, 0.0, 0.0],
            "is_dof_fixed": numpy.full((3,), i == 0, dtype=bool).tolist()
        })

    fixture["rigid_body_problem"]["gravity"] = [0, -9.81, 0]
    fixture["rigid_body_problem"]["coefficient_restitution"] = -1

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
