#!/usr/local/bin/python3
"""
Script to generate a fixture of a chain with N complex links.

Usage: python generate_complex_chainmail_fixture.py N
"""

import json
import pathlib
import sys

import numpy
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

from default_fixture import generate_default_fixture


def create_rotation_matrix(theta: float) -> numpy.ndarray:
    c, s = numpy.cos(theta), numpy.sin(theta)
    return numpy.array(((c, -s), (s, c)))


def create_rectangle(hx: float, hy: float, center: numpy.ndarray, angle: float) -> Polygon:
    points = numpy.array([[hx, hy], [-hx, hy], [-hx, -hy], [hx, -hy]])
    points = points @ create_rotation_matrix(angle).T
    points += center
    return Polygon(points)


def create_link_polygons():
    half_thickness = 1e-2
    width = 6
    height = 5

    head_hx = width / 2 - 4 * half_thickness - 0.5
    foot_hx = width / 4 - 4 * half_thickness
    leg_hy = 2 * height / 7 - half_thickness
    leg_cy = leg_hy + half_thickness
    torso_hx = width / 2
    torso_cy = 2 * leg_cy - half_thickness
    neck_hy = (height - torso_cy - half_thickness) / 2
    neck_cy = height - neck_hy - half_thickness

    return [
        # Head
        create_rectangle(head_hx, half_thickness,
                         numpy.array([width / 2, height - half_thickness]), 0),
        # Neck
        create_rectangle(half_thickness, neck_hy,
                         numpy.array([width / 2, neck_cy]), 0),
        # # Torso
        create_rectangle(torso_hx, half_thickness,
                         numpy.array([torso_hx, torso_cy]), 0),
        # # Left leg
        create_rectangle(half_thickness, leg_hy,
                         numpy.array([half_thickness, leg_cy]), 0),
        # # Left foot
        create_rectangle(foot_hx, half_thickness,
                         numpy.array([foot_hx, half_thickness]), 0),
        # # Right leg
        create_rectangle(half_thickness, leg_hy,
                         numpy.array([width - half_thickness, leg_cy]), 0),
        # # Right foot
        create_rectangle(foot_hx, half_thickness,
                         numpy.array([width - foot_hx, half_thickness]), 0),
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
    R = create_rotation_matrix(theta)

    num_points = vertices.shape[0]
    edges = numpy.hstack([
        numpy.arange(num_points).reshape(-1, 1),
        numpy.roll(numpy.arange(num_points).reshape(-1, 1), -1)]).tolist()

    fixture = generate_default_fixture()
    fixture["timestep_size"] = 1e-2
    fixture["distance_barrier_constraint"]["custom_initial_epsilon"] = 1e-3
    fixture["barrier_solver"]["min_barrier_epsilon"] = 1e-4
    rigid_bodies = fixture["rigid_body_problem"]["rigid_bodies"]
    for i in range(n_links):
        rigid_bodies.append({
            "vertices": vertices.tolist(),
            "edges": edges,
            "polygons": link_polygons,
            "position": (R @ numpy.array([0, -3.5 * i])).tolist(),
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
